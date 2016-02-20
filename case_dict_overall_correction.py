def return_case_df():
    import pandas as pd

    ste_height_correction   = -0.4
    ste_rotation_correction = -1.0
    ste_x_correction        = -2.5

    z00_height_correction   =  0.1
    z00_rotation_correction = -1.0
    z00_x_correction        = -1.0

    z05_height_correction   = -0.2
    z05_rotation_correction = -2.0
    z05_x_correction        = -0.5

    z10_height_correction   =  0.1
    z10_rotation_correction =  0.0
    z10_x_correction        =  0.0

    case_dicts = [ 

        {
            'file':              'STE_a0_p0_U20_z00_tr.hdf5',
            'case_name':         'Straight edge',
            'type':              'surface',
            'x_loc':             "-1,  -5",  # [mm]
            'y_trust_min':       "1, 1", # [mm]
            'Cf':                0.0015,
            'height_correction': ste_height_correction,
            'rotation':          ste_rotation_correction,
            'x_corr':            ste_x_correction,
        },
        {
            'file':              'Sr20R21_a0_p0_U20_z00_tr.hdf5',
            'case_name':         'Serrated, $z/\\lambda = 0$',
            'type':              'surface',
            'x_loc':             "-5, -1, 15, 20, 35, 40", # [mm]
            'y_trust_min':       "0.4, 0.4, 0.4, 0.4, 0.4, 0.4",
            'Cf':                0.0015,
            'height_correction': z00_height_correction,
            'rotation':          z00_rotation_correction,
            'x_corr':            z00_x_correction,
        },
        {
            'file':              'Sr20R21_a0_p0_U20_z05_tr_New.hdf5',
            'case_name':         'Serrated, $z/\\lambda = 0.25$',
            'type':              'surface',
            'x_loc':             "15, 20", # [mm]
            'y_trust_min':       "0.3, 0.3", # [mm]
            'Cf':                0.0015,
            'height_correction': z05_height_correction,
            'rotation':          z05_rotation_correction,
            'x_corr':            z05_x_correction,
        },
        {
            'file':              'Sr20R21_a0_p0_U20_z10_tr.hdf5',
            'case_name':         'Serrated, $z/\\lambda = 0.5$',
            'type':              'surface',
            'x_loc':             "-5, -1", # [mm]
            'y_trust_min':       "1,1", # [mm]
            'Cf':                0.0015,
            'height_correction': z10_height_correction,
            'rotation':          z10_rotation_correction,
            'x_corr':            z10_x_correction,
        },

    ]

    for c in case_dicts:
        if c == case_dicts[0]:
            df = pd.DataFrame( data = c , index = [0])
        else:
            df = df.append( pd.DataFrame( data = c , index = [0]), 
                           ignore_index=True )
    return df
