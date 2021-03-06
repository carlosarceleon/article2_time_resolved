def run_case(case):
    import article2_time_resolved_routines as trr

    trr.raw_data_to_pandas_hdf5(
        case            = case,
        root            = '' ,
        output_root     = '',
        overwrite       = False,
        time_step_limit = 0,
        plot            = False,
        airfoil_normal  = True,
    )

    return 0

import case_dict_overall_correction as case_dict 
case_constants = case_dict.return_case_df()

cases = case_constants.file.values

for c in cases:

    if "STE" in c:
        run_case(c)

