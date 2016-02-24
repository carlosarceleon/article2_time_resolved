

def get_relevant_wall_normal_data_from_pandas_hdf(exceptions = []):
    import article2_time_resolved_routines as trr
    import case_dict_overall_correction as cdoc
    from os.path import join

    root = '/media/carlos/6E34D2CD34D29783/2015-02_SerrationPIV/'\
            +'TR_Data_Location_Calibrated'

    case_dict = cdoc.return_case_df()

    for case in case_dict.iterrows():
        case_file = case[1].file

        skip = False
        if len(exceptions):
            for ex in exceptions:
                if ex in case_file:
                    skip = True

        if not skip:
            x_locs = [float(f) for f in case[1].x_loc.split(',')]

            for x in x_locs:
                if x < 0:
                    cf = case_file.replace('.hdf5','_AirfoilNormal.hdf5')
                else:
                    cf = case_file

                print "  Getting the wall normal data at {0} of \n{1}"\
                        .format( x, cf )
                trr.wall_normal_data_to_reserved_pickles_from_pandas_hdf( 
                    join(root, cf), x , overwrite = False
                )

def get_available_cases_df():
    from os import listdir
    from re import findall
    import pandas as pd

    available_reserved_files = [
        f for f in listdir(root) if f.endswith(".p")
    ]

    cases_df = pd.DataFrame( columns = ['file','x','y'] )

    for res_file in available_reserved_files:
        cases_df = cases_df.append(
            pd.DataFrame( 
                data = {
                    'x' : float(findall( 
                        'x[-0-9][0-9]?.[0-9]', res_file
                    )[0].replace('x','')),

                    'y' : float(findall( 
                        'y[-0-9][0-9]?.?[0-9]?', res_file
                    )[0].replace('y','')),

                    'file' : findall( 
                        '[A-Za-z0-9_]+_x', res_file
                    )[0]\
                    .replace('_x','')\
                    .replace('5','')\
                    .replace('_AirfoilNormal','') ,
                },
                index = [0]
            ),
            ignore_index = True)

    return cases_df

def get_df_cases_from_pickle_names(cases_df):
    from os import listdir
    from os.path import join
    from pandas import read_pickle

    all_pickles = [f for f in listdir(root) if f.endswith('.p')]

    case_time_series_dfs = []
    for case in cases_df.iterrows():
        for pickle in [p for p in all_pickles\
                       if case[1].file in p\
                       and "x{0:.0f}".format(case[1].x) in p\
                       and "y{0:.0f}".format(case[1].y) in p]:

           df = read_pickle( join( root, pickle ) )
           df['case_name'] = case[1].file
           case_time_series_dfs.append( df )

    return case_time_series_dfs


def do_the_time_resolved_analysis():
    import pandas as pd
    from os.path import join
    import article2_data_analysis_routines as dar

    def do_the_frequency_plot(df,plot_name, schematic = ''):
        for y in df.near_y.unique():

            df_y_cases = df[ df.near_y == y ]

            dar.do_the_frequency_analysis( 
                df_y_cases,
                y = y,
                plot_name = plot_name,
                schematic = schematic
            )

    def do_the_Reynolds_quadrant_analysis(df, plot_name):
        for y in TE_cases.near_y.unique():

            df_y_cases = df[ df.near_y == y ]

            dar.do_the_reynolds_stress_quadrant_analysis(
                df_y_cases,
                y = y,
                plot_name = plot_name,
            )

    def do_the_coherence_analysis(df,plot_name,schematic = ''):
        coherence_df = pd.DataFrame()
        for y in TE_cases.near_y.unique():

            partial_coherence_df = dar.do_the_coherence_analysis(
                df[ df.near_y == y ],
            )

            if not partial_coherence_df.empty:
                coherence_df = coherence_df.append( 
                    partial_coherence_df, ignore_index = True 
                )

        dar.plot_coherence_Uc_phi( coherence_df ,
                                 plot_name = plot_name,
                                 schematic = schematic)


    all_cases_pickle = pd.read_pickle( join( root, 'AllPointPickle.p') )

    x0_cases = all_cases_pickle[ 
        all_cases_pickle.near_x == -1 
    ]
    x0_coherence_cases = all_cases_pickle[ 
        (all_cases_pickle.near_x == -3) | (all_cases_pickle.near_x == -1)
    ]

    # TE locations # ###########################################################
    TE_cases = all_cases_pickle[ 
        (all_cases_pickle.near_x == -1) & \
        (all_cases_pickle.case_name == 'STE_a0_p0_U20_z00_tr_AirfoilNormal5')
    ]
    TE_cases = TE_cases.append(
        all_cases_pickle[ 
            (all_cases_pickle.near_x == -1) & \
            (all_cases_pickle.case_name == \
             'Sr20R21_a0_p0_U20_z10_tr_AirfoilNormal5')
        ], ignore_index = True
    )
    TE_cases = TE_cases.append(
        all_cases_pickle[ 
            (all_cases_pickle.near_x == 20) & \
            (all_cases_pickle.case_name == \
             'Sr20R21_a0_p0_U20_z05_tr_New5')
        ], ignore_index = True
    )
    TE_cases = TE_cases.append(
        all_cases_pickle[ 
            (all_cases_pickle.near_x == 40) & \
            (all_cases_pickle.case_name == \
             'Sr20R21_a0_p0_U20_z00_tr5')
        ], ignore_index = True
    )

    # Upwind locations #########################################################
    TE_cases_upwind = all_cases_pickle[ 
        (all_cases_pickle.near_x == -3) & \
        (all_cases_pickle.case_name == 'STE_a0_p0_U20_z00_tr_AirfoilNormal5')
    ]
    TE_cases_upwind = TE_cases_upwind.append(
        all_cases_pickle[ 
            (all_cases_pickle.near_x == -3) & \
            (all_cases_pickle.case_name == \
             'Sr20R21_a0_p0_U20_z10_tr_AirfoilNormal5')
        ], ignore_index = True
    )
    TE_cases_upwind = TE_cases_upwind.append(
        all_cases_pickle[ 
            (all_cases_pickle.near_x == 18) & \
            (all_cases_pickle.case_name == \
             'Sr20R21_a0_p0_U20_z05_tr_New5')
        ], ignore_index = True
    )
    TE_cases_upwind = TE_cases_upwind.append(
        all_cases_pickle[ 
            (all_cases_pickle.near_x == 38) & \
            (all_cases_pickle.case_name == \
             'Sr20R21_a0_p0_U20_z00_tr5')
        ], ignore_index = True
    )

    TE_cases_for_coherence = TE_cases_upwind.append(
        TE_cases, ignore_index = True
    )

    schematic_TE = '/home/carlos/Documents/PhD/Articles/Article_2/'+\
            'Figures/measurement_locations_TE_m2.png'

    schematic_x0 = '/home/carlos/Documents/PhD/Articles/Article_2/'+\
            'Figures/measurement_locations_x0_m2.png'


    do_the_frequency_plot( x0_cases, 'x0', schematic = schematic_x0 )
    do_the_frequency_plot( TE_cases, 'TE',  schematic = schematic_TE)

    #do_the_Reynolds_quadrant_analysis( x0_cases, 'x0' )
    #do_the_Reynolds_quadrant_analysis( TE_cases, 'TE' )

    #do_the_coherence_analysis( x0_coherence_cases , "x0",
    #                         schematic = schematic_x0)
    #do_the_coherence_analysis( TE_cases_for_coherence, 'TE' , 
    #                          schematic = schematic_TE)


from os.path import join
import publish
root = join('/home/carlos/Documents/PhD/Articles/Article_2',
            'Article2_Scripts/time_resolved_scripts/ReservedData/')


get_relevant_wall_normal_data_from_pandas_hdf(exceptions = ['STE'])
get_relevant_wall_normal_data_from_pandas_hdf()
#do_the_time_resolved_analysis()
#publish.publish()
