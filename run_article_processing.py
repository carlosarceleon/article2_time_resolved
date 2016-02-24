##############################################
# Create the HDF5 from in the HD from the
# TECPLOT data
##############################################

def cook_raw_data():
    import os
    import time_data_functions as tdf
    root = "/media/carlos/6E34D2CD34D29783/"\
            +"2015-02_SerrationPIV/TR_Data_NewProcessing/"

    case_folders = [f for f in os.listdir(root) \
                   if os.path.isdir(os.path.join(root,f))\
                   if not os.path.isfile(os.path.join(root,f+".hdf5"))\
                   if not f == 'Parking']

    for cf in case_folders:
       tdf.build_data_hdf5(root,
                           cf,
                           os.path.join(root,cf),
                           overwrite=False
                          )

##############################################
# Grab and pickle all individual points  
##############################################

def build_individual_point_pickles():
    from scitech_manuscript_functions import build_all_cases_popular_lines
    build_all_cases_popular_lines()

##############################################
# Re-pickle all individual points into a
# single pickle
##############################################

def concatenate_individual_point_pickles():
    import scitech_manuscript_functions as proc
    import os
    reload(proc)

    root = './point_data/'

    DFs = proc.all_pickled_coordinates_to_DFs(root)

    for DF in DFs:
        if not os.path.isfile("{0}.p".format(DF.case.unique()[0])):
            DF.to_pickle("{0}.p".format(DF.case.unique()[0]))

def plot_cross_correlation_locations():
    import scitech_manuscript_functions as smf
    reload(smf)
    x_locs    = [0,0.5,0.9]
    y_locs    = [0.1,.3,.6] # Delta normalized
    def make_serration_misalignment_comparison_z05():
        smf.plot_cross_correlation_locations(
            cases     = ['Sr20R21_a0_p0_U20_z05_tr_New',
                         'Sr20R21_a12_p0_U20_z05_tr',
                         'Sr20R21_a12_p6_U20_z05_tr',
                        ],
            case_names = [
                '$\\alpha_g = 0^\\circ,\,\\varphi = 0^\\circ$',
                '$\\alpha_g = 12^\\circ,\,\\varphi = 0^\\circ$',
                '$\\alpha_g = 12^\\circ,\,\\varphi = 6^\\circ$',
            ],
            root      = './WallNormalData/',
            x_locs    = [0,0.5,0.9],
            y_locs    = [0.1,.3,.6], # Delta normalized
            component = 'vx',
            presentation = False,
            plot_name = 'Results/Serration_misalignment_comparison_vx_z05.png'
        )
        smf.plot_cross_correlation_locations(
            cases     = ['Sr20R21_a0_p0_U20_z05_tr_New',
                         'Sr20R21_a12_p0_U20_z05_tr',
                         'Sr20R21_a12_p6_U20_z05_tr',
                         
                        ],
            case_names = [
                '$\\alpha_g = 0^\\circ,\,\\varphi = 0^\\circ$',
                '$\\alpha_g = 12^\\circ,\,\\varphi = 0^\\circ$',
                '$\\alpha_g = 12^\\circ,\,\\varphi = 6^\\circ$',
            ],
            root      = './WallNormalData/',
            x_locs    = x_locs,
            y_locs    = y_locs,
            component = 'vy',
            presentation = False,
            plot_name = 'Results/Serration_misalignment_comparison_vy_z05.png'
        )
    def make_serration_misalignment_comparison_z10():
        smf.plot_cross_correlation_locations(
            cases     = ['Sr20R21_a0_p0_U20_z10_tr',
                         'Sr20R21_a12_p0_U20_z10_tr_old',
                         'Sr20R21_a12_p6_U20_z10_tr',
                        ],
            case_names = [
                '$\\alpha_g = 0^\\circ,\,\\varphi = 0^\\circ$',
                '$\\alpha_g = 12^\\circ,\,\\varphi = 0^\\circ$',
                '$\\alpha_g = 12^\\circ,\,\\varphi = 6^\\circ$',
            ],
            root      = './WallNormalData/',
            x_locs    = x_locs,
            y_locs    = y_locs,
            component = 'vx',
            presentation = False,
            plot_name = 'Results/Serration_misalignment_comparison_vx_z10.png'
        )
        smf.plot_cross_correlation_locations(
            cases     = ['Sr20R21_a0_p0_U20_z10_tr',
                         'Sr20R21_a12_p0_U20_z10_tr_old',
                         'Sr20R21_a12_p6_U20_z10_tr',
                        ],
            case_names = [
                '$\\alpha_g = 0^\\circ,\,\\varphi = 0^\\circ$',
                '$\\alpha_g = 12^\\circ,\,\\varphi = 0^\\circ$',
                '$\\alpha_g = 12^\\circ,\,\\varphi = 6^\\circ$',
            ],
            root      = './WallNormalData/',
            x_locs    = x_locs,
            y_locs    = y_locs,
            component = 'vy',
            presentation = False,
            plot_name = 'Results/Serration_misalignment_comparison_vy_z10.png'
        )
    def make_serration_spanwise_comparison_no_miss():
        smf.plot_cross_correlation_locations(
            cases     = [
                         'STE_a0_p0_U20_z00_tr',
                         'Sr20R21_a0_p0_U20_z00_tr',
                         'Sr20R21_a0_p0_U20_z05_tr_New',
                         'Sr20R21_a0_p0_U20_z10_tr',
                        ],
            case_names= [
                'Straight trailing edge',
                'Serrated, $z/\\lambda=0$',
                'Serrated, $z/\\lambda=0.25$',
                'Serrated, $z/\\lambda=0.5$',
            ],
            root      = './WallNormalData/',
            x_locs    = x_locs,
            y_locs    = y_locs,
            component = 'vx',
            presentation = False,
            plot_name = 'Results/Serration_no_misalignment_comparison_vx.png'
        )
        smf.plot_cross_correlation_locations(
            cases     = [
                         'STE_a0_p0_U20_z00_tr',
                         'Sr20R21_a0_p0_U20_z00_tr',
                         'Sr20R21_a0_p0_U20_z05_tr_New',
                         'Sr20R21_a0_p0_U20_z10_tr',
                        ],
            case_names= [
                'Straight trailing edge',
                'Serrated, $z/\\lambda=0$',
                'Serrated, $z/\\lambda=0.25$',
                'Serrated, $z/\\lambda=0.5$',
            ],
            root      = './WallNormalData/',
            x_locs    = x_locs,
            y_locs    = y_locs,
            component = 'vy',
            presentation = False,
            plot_name = 'Results/Serration_no_misalignment_comparison_vy.png'
        )

    def make_new_z05():
        smf.plot_cross_correlation_locations(
            cases     = [
                         'Sr20R21_a0_p0_U20_z05_tr',
                         'Sr20R21_a0_p0_U20_z05_tr_New',
                        ],
            root      = './WallNormalData/',
            x_locs    = x_locs,
            y_locs    = y_locs,
            component = 'vx',
            presentation = False,
            plot_name = 'Results/New_z05_vx.png'
        )
        smf.plot_cross_correlation_locations(
            cases     = [
                         'Sr20R21_a0_p0_U20_z05_tr',
                         'Sr20R21_a0_p0_U20_z05_tr_New',
                        ],
            root      = './WallNormalData/',
            x_locs    = x_locs,
            y_locs    = y_locs,
            component = 'vy',
            presentation = False,
            plot_name = 'Results/New_z05_vy.png'
        )
    def test():
        smf.plot_cross_correlation_locations(
            cases     = [
                         'Sr20R21_a0_p0_U20_z00_tr',
                        ],
            root      = './WallNormalData/',
            x_locs    = x_locs,
            y_locs    = y_locs,
            component = 'vx',
            presentation = False,
            plot_name = 'Results/test.png',
            test      = True
        )

    make_serration_spanwise_comparison_no_miss()
    make_serration_misalignment_comparison_z05()
    make_serration_misalignment_comparison_z10()
    #make_new_z05()
    #test()


##############################################
# Line PSD plot of interesting cases, 
# locations
##############################################
def plot_select_frequency_locations():
    import scitech_manuscript_functions as smf
    reload(smf)
    def make_serration_misalignment_comparison_z10():
        smf.plot_select_frequency_locations(
            cases     = [
                'Sr20R21_a0_p0_U20_z10_tr',
                'Sr20R21_a12_p0_U20_z10_tr_old',
                'Sr20R21_a12_p6_U20_z10_tr',
            ],
            case_names = [
                '$\\alpha_g = 0^\\circ,\\,\\varphi = 0^\\circ$',
                '$\\alpha_g = 12^\\circ,\\,\\varphi = 0^\\circ$',
                '$\\alpha_g = 12^\\circ,\\,\\varphi = 6^\\circ$',
            ],
            root      = './WallNormalData/',
            x_locs    = [0   , 0.5 , 1.] ,
            y_locs    = [0.1 , .3  , .6] , # Delta normalized
            component = 'vx',
            presentation = False,
            plot_name = 'SpectrumResults/Serration_misalignment_comparison_vx_z10.png'
        )
        smf.plot_select_frequency_locations(
            cases     = ['Sr20R21_a12_p0_U20_z10_tr_old',
                         'Sr20R21_a12_p6_U20_z10_tr',
                         'Sr20R21_a0_p0_U20_z10_tr',
                        ],
            case_names = [
                '$\\alpha_g = 0^\\circ,\\,\\varphi = 0^\\circ$',
                '$\\alpha_g = 12^\\circ,\\,\\varphi = 0^\\circ$',
                '$\\alpha_g = 12^\\circ,\\,\\varphi = 6^\\circ$',
            ],
            root      = './WallNormalData/',
            x_locs    = [0,0.5,1.],
            y_locs    = [0.1,.3,.6], # Delta normalized
            component = 'vy',
            presentation = False,
            plot_name = 'SpectrumResults/Serration_misalignment_comparison_vy_z10.png'
        )

    def make_serration_misalignment_comparison_z05():
        smf.plot_select_frequency_locations(
            cases     = ['Sr20R21_a0_p0_U20_z05_tr_New',
                         'Sr20R21_a12_p0_U20_z05_tr',
                         'Sr20R21_a12_p6_U20_z05_tr',
                         
                        ],
            case_names = [
                '$\\alpha_g = 0^\\circ,\,\\varphi = 0^\\circ$',
                '$\\alpha_g = 12^\\circ,\,\\varphi = 0^\\circ$',
                '$\\alpha_g = 12^\\circ,\,\\varphi = 6^\\circ$',
            ],
            root      = './WallNormalData/',
            x_locs    = [0   , 0.5 , 1.] ,
            y_locs    = [0.1 , .3  , .6] , # Delta normalized
            component = 'vx',
            presentation = False,
            plot_name = 'SpectrumResults/Serration_misalignment_comparison_vx_z05.png'
        )
        smf.plot_select_frequency_locations(
            cases     = [
                'Sr20R21_a0_p0_U20_z05_tr_New',
                'Sr20R21_a12_p0_U20_z05_tr',
                'Sr20R21_a12_p6_U20_z05_tr',
                        ],
            case_names = [
                '$\\alpha_g = 0^\\circ,\,\\varphi = 0^\\circ$',
                '$\\alpha_g = 12^\\circ,\,\\varphi = 0^\\circ$',
                '$\\alpha_g = 12^\\circ,\,\\varphi = 6^\\circ$',
            ],
            root      = './WallNormalData/',
            x_locs    = [0,0.5,1.],
            y_locs    = [0.1,.3,.6], # Delta normalized
            component = 'vy',
            presentation = False,
            plot_name = 'SpectrumResults/Serration_misalignment_comparison_vy_z05.png'
        )
    def make_serration_spanwise_comparison_no_miss():
        smf.plot_select_frequency_locations(
            cases     = [
                         'STE_a0_p0_U20_z00_tr',
                         'Sr20R21_a0_p0_U20_z00_tr',
                         'Sr20R21_a0_p0_U20_z05_tr_New',
                         'Sr20R21_a0_p0_U20_z10_tr',
                        ],
            case_names= [
                'Straight trailing edge',
                'Serrated, $z/\\lambda=0$',
                'Serrated, $z/\\lambda=0.25$',
                'Serrated, $z/\\lambda=0.5$',
            ],
            root      = './WallNormalData/',
            x_locs    = [0,0.5,1.],
            y_locs    = [0.1,.3,.6], # Delta normalized
            component = 'vx',
            presentation = False,
            plot_name = 'SpectrumResults/Serration_no_misalignment_comparison_vx.png'
        )
        smf.plot_select_frequency_locations(
            cases     = [
                         'STE_a0_p0_U20_z00_tr',
                         'Sr20R21_a0_p0_U20_z00_tr',
                         'Sr20R21_a0_p0_U20_z05_tr_New',
                         'Sr20R21_a0_p0_U20_z10_tr',
                        ],
            case_names= [
                'Straight trailing edge',
                'Serrated, $z/\\lambda=0$',
                'Serrated, $z/\\lambda=0.25$',
                'Serrated, $z/\\lambda=0.5$',
            ],
            root      = './WallNormalData/',
            x_locs    = [0,0.5,1.],
            y_locs    = [0.1,.3,.6], # Delta normalized
            component = 'vy',
            presentation = False,
            plot_name = 'SpectrumResults/Serration_no_misalignment_comparison_vy.png'
        )

    def make_serration_ste_comparison():
        smf.plot_select_frequency_locations(
            cases     = [
                         'Sr20R21_a0_p0_U20_z00_tr',
                         'STE_a0_p0_U20_z00_tr',
                        ],
            case_names = [
                'Straight',
                'Serrated',
            ],
            root      = './WallNormalData/',
            x_locs    = [0,0.5,1.],
            y_locs    = [0.1,.3,.6], # Delta normalized
            component = 'vx',
            presentation = False,
            plot_name = 'SpectrumResults/Serration_STE_comparison_vx.png'
        )
        smf.plot_select_frequency_locations(
            cases     = [
                         'Sr20R21_a0_p0_U20_z00_tr',
                         'STE_a0_p0_U20_z00_tr',
                        ],
            case_names = [
                'Straight',
                'Serrated',
            ],
            root      = './WallNormalData/',
            x_locs    = [0,0.5,1.],
            y_locs    = [0.1,.3,.6], # Delta normalized
            component = 'vy',
            presentation = False,
            plot_name = 'SpectrumResults/Serration_STE_comparison_vy.png'
        )
    def make_new_z05():
        smf.plot_select_frequency_locations(
            cases     = [
                         'Sr20R21_a0_p0_U20_z05_tr',
                         'Sr20R21_a0_p0_U20_z05_tr_New',
                        ],
            root      = './WallNormalData/',
            x_locs    = [0,0.5,1.],
            y_locs    = [0.1,.3,.6], # Delta normalized
            component = 'vx',
            presentation = False,
            plot_name = 'SpectrumResults/New_z05_vx.png'
        )
        smf.plot_select_frequency_locations(
            cases     = [
                         'Sr20R21_a0_p0_U20_z05_tr',
                         'Sr20R21_a0_p0_U20_z05_tr_New',
                        ],
            root      = './WallNormalData/',
            x_locs    = [0,0.5,1.],
            y_locs    = [0.1,.3,.6], # Delta normalized
            component = 'vy',
            presentation = False,
            plot_name = 'SpectrumResults/New_z05_vy.png'
        )
    make_serration_spanwise_comparison_no_miss()
    make_serration_misalignment_comparison_z10()
    make_serration_ste_comparison()
    make_serration_misalignment_comparison_z05()
    #make_new_z05()

def make_freq_maps():

    import scitech_manuscript_functions as tr_func
    import pandas as pd
    import os
    reload(tr_func)

    root = '/home/carlos/Documents/PhD/Articles/Article_2'+\
            '/Article2_Scripts/time_resolved_scripts/point_data'

    pickled_cases = [f for f in os.listdir('.') \
                     if f.endswith('.p')]

    for pc in pickled_cases:
        if not os.path.isfile(
            "./Presentation_FrequencyMap"+pc.replace('.p','.png')
        ):
            df = pd.read_pickle(os.path.join('.',pc))

            tr_func.make_wall_normal_frequency_map(
                case_name    = pc.replace('.p',''),
                root         = root,
                df           = df,
                plot_name    = "./FrequencyMap"+pc.replace('.p',''),
                presentation = True
            )

def check_timestep():
    import time_data_functions as tdf
    import os
    reload(tdf)

    device   = 'STE'
    alpha    = '0'
    phi      = '0'
    z        = '00'
    t        = [0]
    variable = 'vy'

    root = "/media/carlos/6E34D2CD34D29783/"\
            +"2015-02_SerrationPIV/TR_Data_NewProcessing/"

    case_name = "{0}_a{1}_p{2}_U20_z{3}_tr".format(device,alpha,phi,z)
    hdf5_file = os.path.join(root,case_name+".hdf5")
    plt_name = case_name+"_t{0:.2f}".format(t[0])+".png"


    tdf.plot_surface(
        case      = case_name,
        hdf5_file = hdf5_file,
        plt_name  = plt_name,
        time      = t,
        variable  = variable
    )

def extract_wall_normal_lines():
    from subprocess import call

    call(['python','run_article_processing.py'])

def new_rotated_hdf5_file():
    import os
    import time_data_functions as tdf

    root = "/media/carlos/6E34D2CD34D29783/"\
            +"2015-02_SerrationPIV/TR_Data_NewProcessing/"

    cases = [f.replace('.hdf5','') for f in os.listdir(root) \
             if f.endswith('.hdf5')]
             
    for case_name in [cases[0]]:
        tdf.return_rotated_hdf5_in_df(
                hdf5_file     = os.path.join(root,case_name+".hdf5") ,
                case          = case_name,
                output_name   = case_name+'_Rotated.hdf5',
                overwrite     = False,
            )

def extract_wall_normal_lines_one_by_one():

    import os
    import time_data_functions as tdf

    root = "/media/carlos/6E34D2CD34D29783/"\
            +"2015-02_SerrationPIV/TR_Data_NewProcessing/"

    cases = [f.replace('.hdf5','') for f in os.listdir(root) \
             if f.endswith('.hdf5')]
             
    for case_name in cases:
        tdf.get_time_resolved_wall_normal_line(
            hdf5_file     = os.path.join(root,case_name+".hdf5") ,
            case          = case_name,
            output_hdf    = os.path.join(
                'WallNormalData' , case_name+'_WallNormalData.hdf5'
            ),
            plot          = False,
            overwrite     = False,
            x_locs        = [
                -0.10 , -0.05, -0.025 , 0.025, 0.05  , 0.10 ,
                0.20  , 0.30  ,
                0.40  , 0.5   , 0.60 ,
                0.70  , 0.80  ,
                0.90  , 1.0   , 1.10 ,
            ],
        )


#cook_raw_data()
#extract_wall_normal_lines_one_by_one()
#make_freq_maps()
#plot_select_frequency_locations()
plot_cross_correlation_locations()
#check_timestep()

