

def get_relevant_wall_normal_data_from_pandas_hdf():
    import article2_time_resolved_routines as trr
    import case_dict_overall_correction as cdoc
    from os.path import join

    root = '/media/carlos/6E34D2CD34D29783/2015-02_SerrationPIV/'\
            +'TR_Data_Location_Calibrated'

    case_dict = cdoc.return_case_df()

    for case in case_dict.iterrows():
        case_file = case[1].file
        x_locs = [float(f) for f in case[1].x_loc.split(',')]

        for x in x_locs:
            if x < 0:
                cf = case_file.replace('.hdf5','_AirfoilNormal.hdf5')
            else:
                cf = case_file

            print "  Getting the wall normal data at {0} of \n{1}"\
                    .format( x, cf )
            trr.wall_normal_data_to_reserved_pickles_from_pandas_hdf( 
                join(root, cf), x 
            )

get_relevant_wall_normal_data_from_pandas_hdf()
