def run_case(case_name):
    import time_data_functions as tdf
    import os
    root = "/media/carlos/6E34D2CD34D29783/"\
            +"2015-02_SerrationPIV/TR_Data_NewProcessing/"

    print "Submitting {0}".format(case_name)
    tdf.get_time_resolved_wall_normal_line(
        hdf5_file     = os.path.join(root,case_name+".hdf5") ,
        case          = case_name,
        output_hdf    = os.path.join(
            'WallNormalData' , case_name+'_WallNormalData.hdf5'
        ),
        plot          = False,
        overwrite     = False,
        x_locs        = [
            -0.10 , -0.05 , 0.0  , 0.10 ,
            0.20  , 0.30  ,
            0.40  , 0.5   , 0.60 ,
            0.70  , 0.80  ,
            0.90  , 1.0   , 1.10 ,
        ],
    )
    return 0

from multiprocessing import Pool
import os

p = Pool(4)

root = "/media/carlos/6E34D2CD34D29783/"\
        +"2015-02_SerrationPIV/TR_Data_NewProcessing/"

cases = [f.replace('.hdf5','') for f in os.listdir(root) \
         if f.endswith('.hdf5')
        ]
         
p.map(run_case,cases)

