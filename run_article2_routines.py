

def get_relevant_wall_normal_data_from_pandas_hdf(exceptions = []):
    import article2_time_resolved_routines as trr
    import case_dict_overall_correction as cdoc
    from os.path import join
    from numpy import arange

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
            #x_locs = [float(f) for f in case[1].x_loc.split(',')]
            if 'z00' in case_file and 'Sr' in case_file:
                x_locs = arange( 36, 40.2, 0.2 )
            elif 'z05' in case_file and 'Sr' in case_file:
                x_locs = arange( 16, 20.2, 0.2 )
            elif 'z10' in case_file or 'STE' in case_file:
                x_locs = arange( -5, -0.8, 0.2 )

            trr.wall_normal_data_to_reserved_pickles_from_pandas_hdf( 
                join(root, case_file), x_locs , overwrite = False, append = True
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

#def correct_heights(df):
#    from numpy import round as np_round
#    
#    def round_of_rating(number):
#        """Round a number to the closest half integer.
#        """
#        return round(number * 20) / 2 / 10.
#
#    height_correction = {
#        'Sr20R21_a0_p0_U20_z00_tr':     -1.0,
#        'Sr20R21_a0_p0_U20_z05_tr_New':  0.0,
#        'Sr20R21_a0_p0_U20_z10_tr':      1.0,
#    }
#
#    for key in height_correction.keys():
#
#        df.loc[ df.case_name == key , 'y'] = \
#                df[ df.case_name == key ].y + height_correction[key]
#
#        df.loc[ df.case_name == key , 'near_y'] = \
#                df[ df.case_name == key ].near_y \
#                + height_correction[key]
#
#        df.loc[ df.case_name == key , 'near_y_delta'] = \
#                df[ df.case_name == key ].near_y \
#                / df[ df.case_name == key ].delta_99
#
#    df.near_y_delta = np_round(df.near_y_delta,1)
#    df.near_y       = np_round(df.near_y      ,1)
#
#    print sorted(df.near_y_delta.unique())
#    return df


def do_the_time_resolved_analysis():
    import pandas as pd
    from os.path import join
    import article2_data_analysis_routines as dar

    def do_the_frequency_plot(df,plot_name, schematic = ''):
        for y in df.near_y_delta.unique():

            df_y_cases = df[ df.near_y_delta == y ]

            dar.do_the_frequency_analysis( 
                df_y_cases,
                y = y,
                plot_name = plot_name,
                schematic = schematic
            )

    def do_the_Reynolds_quadrant_analysis(df, plot_name):
        for y in df.near_y_delta.unique():

            df_y_cases = df[ df.near_y_delta == y ]

            dar.do_the_reynolds_stress_quadrant_analysis(
                df_y_cases,
                y_delta = y,
                plot_name = plot_name,
            )

    def do_the_coherence_analysis(df_upstream, df_downstream
                                  ,plot_name,schematic = ''):
        coherence_df = pd.DataFrame()
        for y_up, y_down in zip(
            sorted(df_upstream.near_y_delta.unique()),
            sorted(df_downstream.near_y_delta.unique())
        ):

            partial_coherence_df = dar.do_the_coherence_analysis(
                df_upstream[ df_upstream.near_y_delta == y_up ],
                df_downstream[ df_downstream.near_y_delta == y_down ],
            )

            if not partial_coherence_df.empty:
                coherence_df = coherence_df.append( 
                    partial_coherence_df, ignore_index = True 
                )

        dar.plot_coherence_Uc_phi( coherence_df ,
                                 plot_name = plot_name,
                                 schematic = schematic)


    all_cases_pickle = pd.read_pickle( join( root, 'AllPointPickle.p' ) )

    x0_cases = all_cases_pickle[ 
        all_cases_pickle.near_x == -1 
    ]
    x0_coherence_cases = all_cases_pickle[ 
        (all_cases_pickle.near_x == -3) | (all_cases_pickle.near_x == -1)
    ]

    x0_coherence_cases = x0_coherence_cases[
        x0_coherence_cases.case_name != "STE_a0_p0_U20_z00_tr"
    ]

    # TE locations # ###########################################################
    #TE_cases = all_cases_pickle[ 
    #    (all_cases_pickle.near_x == -1) & \
    #    (all_cases_pickle.case_name == 'STE_a0_p0_U20_z00_tr')
    #]
    TE_cases = \
        all_cases_pickle[ 
            #(all_cases_pickle.near_x == -1) & \
            (all_cases_pickle.case_name == \
             'Sr20R21_a0_p0_U20_z10_tr')
        ]
    TE_cases = TE_cases.append(
        all_cases_pickle[ 
            #(all_cases_pickle.near_x == 20) & \
            (all_cases_pickle.case_name == \
             'Sr20R21_a0_p0_U20_z05_tr_New')
        ], ignore_index = True
    )
    TE_cases = TE_cases.append(
        all_cases_pickle[ 
            #(all_cases_pickle.near_x == 40) & \
            (all_cases_pickle.case_name == \
             'Sr20R21_a0_p0_U20_z00_tr')
        ], ignore_index = True
    )
    TE_cases = TE_cases.append(
        all_cases_pickle[ 
            #(all_cases_pickle.near_x == -1) & \
            (all_cases_pickle.case_name == \
             'STE_a0_p0_U20_z00_tr')
        ], ignore_index = True
    )


    # upstream locations ######################################################
    #TE_cases_upstream = all_cases_pickle[ 
    #    (all_cases_pickle.near_x == -3) & \
    #    (all_cases_pickle.case_name == 'STE_a0_p0_U20_z00_tr')
    #]
    #up_shift = -4
    up_shift = -2
    TE_cases_upstream = \
        all_cases_pickle[ 
            (all_cases_pickle.near_x == -1 + up_shift) & \
            (all_cases_pickle.case_name == \
             'Sr20R21_a0_p0_U20_z10_tr')
        ]
    TE_cases_upstream = TE_cases_upstream.append(
        all_cases_pickle[ 
            (all_cases_pickle.near_x == 20 + up_shift) & \
            (all_cases_pickle.case_name == \
             'Sr20R21_a0_p0_U20_z05_tr_New')
        ], ignore_index = True
    )
    TE_cases_upstream = TE_cases_upstream.append(
        all_cases_pickle[ 
            (all_cases_pickle.near_x == 40 + up_shift) & \
            (all_cases_pickle.case_name == \
             'Sr20R21_a0_p0_U20_z00_tr')
        ], ignore_index = True
    )
    TE_cases_upstream = TE_cases_upstream.append(
        all_cases_pickle[ 
            (all_cases_pickle.near_x == -1 + up_shift) & \
            (all_cases_pickle.case_name == \
             'STE_a0_p0_U20_z00_tr')
        ], ignore_index = True
    )

    #TE_cases          = correct_heights( TE_cases )
    #TE_cases_upstream = correct_heights( TE_cases_upstream )


    schematic_TE = '/home/carlos/Documents/PhD/Articles/Article_2/'+\
            'Figures/measurement_locations_TE_m2_with_edge_normal.png'

    schematic_x0 = '/home/carlos/Documents/PhD/Articles/Article_2/'+\
            'Figures/measurement_locations_x0_m2.png'


    schematic_TE = '/home/carlos/Documents/PhD/Articles/Article_2/'+\
            'Figures/measurement_locations_TE_m2_noSTE.png'

    #TE_cases = TE_cases[ TE_cases.case_name != 'STE_a0_p0_U20_z00_tr' ]

    #do_the_Reynolds_quadrant_analysis( TE_cases, 'TE' )

    #do_the_coherence_analysis( TE_cases_upstream, TE_cases, 'TE' , 
    #                          schematic = schematic_TE)
    #dar.plot_mean_and_std( TE_cases )

    #do_the_frequency_plot( x0_cases, 'x0', schematic = schematic_x0 )
    #do_the_Reynolds_quadrant_analysis( x0_cases, 'x0' )
    #do_the_coherence_analysis( x0_coherence_cases , "x0",
    #                         schematic = schematic_x0)

def correlation_coherence_and_length_scale_analysis():
    import article2_data_analysis_routines as dar

    def do_the_vertical_coherence_analysis( df ):

        dar.plot_vertical_coherence( df )

    def do_the_streamwise_coherence_analysis( pickle ):

        dar.plot_streamwise_f_coherence( pickle )



    #do_the_frequency_plot( TE_cases, 'TE',  schematic = schematic_TE)
    #do_the_streamwise_coherence_analysis( 
    #    '/home/carlos/Documents/PhD/Articles/Article_3/Scripts/' + \
    #    'time_resolved/ReservedData/Sr20R21_a0_p0_U20_z05_tr.p',
    #    'TE' )
    #do_the_streamwise_coherence_analysis( 
    #    '/home/carlos/Documents/PhD/Articles/Article_3/Scripts/' + \
    #    'time_resolved/ReservedData/Sr20R21_a0_p0_U20_z00_tr.p',
    #    'TE' )
    #do_the_streamwise_coherence_analysis( 
    #    '/home/carlos/Documents/PhD/Articles/Article_3/Scripts/' + \
    #    'time_resolved/ReservedData/Sr20R21_a0_p0_U20_z10_tr.p',
    #    'TE' )
    #do_the_streamwise_coherence_analysis( 
    #    '/home/carlos/Documents/PhD/Articles/Article_3/Scripts/' + \
    #    'time_resolved/ReservedData/STE_a0_p0_U20_z00_tr.p',
    #    'TE' )

    #do_the_vertical_coherence_analysis(
    #    [
    #        '/home/carlos/Documents/PhD/Articles/Article_3/Scripts/' + \
    #        'time_resolved/ReservedData/STE_a0_p0_U20_z00_tr.p',

    #        '/home/carlos/Documents/PhD/Articles/Article_3/Scripts/' + \
    #        'time_resolved/ReservedData/Sr20R21_a0_p0_U20_z00_tr.p',

    #        '/home/carlos/Documents/PhD/Articles/Article_3/Scripts/' + \
    #        'time_resolved/ReservedData/Sr20R21_a0_p0_U20_z05_tr.p',

    #        '/home/carlos/Documents/PhD/Articles/Article_3/Scripts/' + \
    #        'time_resolved/ReservedData/Sr20R21_a0_p0_U20_z10_tr.p',

    #    ]
    #    )

    dar.get_vertical_length_scale()


from os.path import join
import publish
root = join('/home/carlos/Documents/PhD/Articles/Article_2',
            'Article2_Scripts/time_resolved_scripts/LineReservedData/')


#get_relevant_wall_normal_data_from_pandas_hdf(exceptions = ['STE'])
#get_relevant_wall_normal_data_from_pandas_hdf()
#get_relevant_wall_normal_data_from_pandas_hdf(exceptions = ['z05','STE','z00'])
#do_the_time_resolved_analysis()
correlation_coherence_and_length_scale_analysis()
#publish.publish()
