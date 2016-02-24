def run_case(case):
    import article2_time_resolved_routines as trr
    import run_article2_routines as routines

    trr.raw_data_to_pandas_hdf5(
        case            = case,
        root            = '' ,
        output_root     = '',
        overwrite       = False,
        time_step_limit = 0,
        plot            = False,
        airfoil_normal  = False,
    )

    trr.raw_data_to_pandas_hdf5(
        case            = case,
        root            = '' ,
        output_root     = '',
        overwrite       = False,
        time_step_limit = 0,
        plot            = False,
        airfoil_normal  = True,
    )

    try:
        routines.get_relevant_wall_normal_data_from_pandas_hdf()
    except:
        pass

    return 0

from multiprocessing import Pool
import case_dict_overall_correction as case_dict 
case_constants = case_dict.return_case_df()

cases = case_constants.file.values

p = Pool(4)

p.map(run_case,cases)

