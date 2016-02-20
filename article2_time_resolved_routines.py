def get_vorticity(df):
    from numpy import shape,zeros
    from sys import exit

    if "vorticity_xy" in df.columns:
        # Do nothing and return the same DF
        return df

    # Get shape of 2D plane
    nx = len(df['x'].unique())
    ny = len(df['y'].unique())
    Ux = df['u'].values.reshape((ny,nx))
    Uy = df['v'].values.reshape((ny,nx))
    ax = df['x'].values.reshape((ny,nx))/1000. # [mm] -> [m]
    ay = df['y'].values.reshape((ny,nx))/1000. # [mm] -> [m]

    i,j = shape(Ux)

    # Sanity check:
    if i != shape(Uy)[0] or i != shape(ax)[0] or i != shape(ay)[0]:
        exit("   The shape of the matrices while getting the "+\
             "vorticity is not the same!")

    duy_dax = zeros((i,j))
    dux_day = zeros((i,j))
    for ii in range(1,i-1):
        for jj in range(1,j-1):
            duy_dax[ii,jj] = (Uy[ii,jj+1]-Uy[ii,jj-1])\
                    /(ax[ii,jj+1]-ax[ii,jj-1])
    for ii in range(1,i-1):
        for jj in range(1,j-1):
            dux_day[ii,jj] = (Ux[ii+1,jj]-Ux[ii-1,jj])\
                    /(ay[ii+1,jj]-ay[ii-1,jj])

    vorticity = duy_dax - dux_day

    df['vorticity_xy'] = vorticity.ravel()

    return df

def read_raw_hdf5_case_and_write_pandas_hdf5(
    hdf5_file,
    root                  = ''    ,
    output_file           = ''    ,
    serration_angle       = 0     ,
    angle_correction      = 0     ,
    height_correction     = 0     ,
    streamwise_correction = 0     ,
    overwrite             = False ,
    time_step_limit       = 0     ,
    plot                  = False ,
    airfoil_normal        = False,
):
    """ Gets an HDF5 file and rotates it into a new one

    Input:
        hdf5_file,
        case: case name
        output_name
        overwrite
    """

    #######################################################
    #######################################################
    # IMPORTANT
    #
    # The coordinates coming from the HDF5 file are the
    # vertical freestream coordinates of DaVis.
    #
    # The coordinates used for the local variables are
    # already put to the left-to-right freestream 
    # coordinates
    #
    #######################################################
    #######################################################

    from progressbar import ProgressBar,Percentage,Bar,ETA,SimpleProgress
    import h5py
    import numpy as np
    import pandas as pd
    from os.path import isfile,join

    write_frequency = 150

    case = hdf5_file.replace('.hdf5','')

    # File related things ######################################################
    if not output_file:
        output_file = case+".hdf5"

    if airfoil_normal:
        output_file = output_file+"_AirfoilNormal"

    if not output_file.endswith('.hdf5'):
        output_file = output_file.replace(".hdf5","")+".hdf5"

    if isfile(output_file) and not overwrite:
        print "  Exiting; file exists:\n{0}".format(output_file)
        return 0
    # ##########################################################################

    h5 = h5py.File(join(root,hdf5_file),'r')

    # Read the available times #################################################
    available_times = sorted([int(f[0]) for f in \
                              h5['{0}'.format(case)].iteritems()\
                              if not 'mask' in f and not 'x' in f and not 'y'\
                              in f])
    # ##########################################################################

    if time_step_limit:
        available_times = available_times[:time_step_limit]

    progress = ProgressBar(
         widgets=[
             Bar(),' ',
             Percentage(),' ',
             ETA(), ' (time step ',
             SimpleProgress(),')'], 
         maxval=len(available_times)
         ).start()

    t_x_cnt = 0
    cnt     = 0

    hdf = pd.HDFStore(output_file)

    df_dump = pd.DataFrame( columns = ['x','y','u','v','w','time_step'] )

    rotation_angle = serration_angle + angle_correction
    if airfoil_normal:
        rotation_angle = rotation_angle - 11.4

    for ti in available_times:
        df = pd.DataFrame( data = {
            'x' :  np.array(h5["{0}/y".format(case)].value),
            'y' : -np.array(h5["{0}/x".format(case)].value),
            'u' :  np.array(h5["{0}/{1}/{2}".format(case,ti,'Vy')].value),
            'v' : -np.array(h5["{0}/{1}/{2}".format(case,ti,'Vx')].value),
            'w' :  np.array(h5["{0}/{1}/{2}".format(case,ti,'Vz')].value),
        })

        df[ 'time_step' ] = ti

        df = correct_flow_plane_df(
            df,
            rotation_angle        = rotation_angle,
            height_correction     = height_correction,
            streamwise_correction = streamwise_correction,
        )

        progress.update(ti)

        df_dump = df_dump.append(df,ignore_index=True)

        if cnt == write_frequency:

            if plot:
                show_surface_from_df(
                    df[df.time_step == ti], 
                    'u'
                )

            if t_x_cnt == cnt:
                hdf.put(
                    case, 
                    df_dump.convert_objects(), 
                    format='table', 
                    data_columns=True
                )

            else:
                hdf.append(
                    case , 
                    df_dump.convert_objects(), 
                    format='table', 
                    data_columns=True
                )

            df_dump = pd.DataFrame( 
                columns = ['x','y','u','v','w','time_step'] 
            )
            cnt = 0

        if ti == available_times[-1]:
            hdf.append(
                case , 
                df_dump.convert_objects(), 
                format='table', 
                data_columns=True
            )

        t_x_cnt += 1
        cnt     += 1

    hdf.close()
    h5.close()

    progress.finish()

def show_surface_from_df(df, variable='u', points = [], mask = []):
    import matplotlib.pyplot as plt
    from numpy import meshgrid,linspace

    if len(df.x.unique())>1000:
        df = regrid_df(df)

    pivoted_df = df.pivot_table(variable, 'x', 'y', fill_value=0)

    y = pivoted_df.columns
    x = pivoted_df.index

    X,Y = meshgrid( x,y )
    Z   = pivoted_df.values.T

    fig,axes = plt.subplots(1,1)

    if variable == 'flow_angle':
        df.flow_angle.fillna(90)
        levels = linspace(-15,15)
        cf = axes.contourf(X,Y,Z,levels=levels)

    else:
        cf = axes.contourf(X,Y,Z)

    if len(mask):
        axes.plot(mask[0,:],mask[1,:])

    if len(points):
        axes.scatter(points[0],points[1])

    plt.colorbar(cf)
    axes.set_aspect('equal')

    plt.show()
    plt.close(fig)

def get_case_details_from_filename(case_name):
    from re import findall
    device  = case_name.split("_")[0]
    phi   = findall("p-?[0-9]"       , case_name )[0].replace("p"   , "")
    alpha = findall("a-?[0-9][0-9]?" , case_name )[0].replace("a" , "")
    U     = findall("U[0-9]+"        , case_name )[0].replace("U"     , "")
    try:
        loc   = findall("z[0-9][0-9]?"   , case_name )[0].replace("z"     , "")
    except IndexError:
        loc = 0
    reprocessed = False
    if "_Reprocessed" in case_name:
       reprocessed = True

    return device,phi,alpha,U,loc,reprocessed

def find_nearest(to_point,from_array):
   """ Finds the nearest available value in a array to a given value

   Inputs:
      to_point: value to find the nearest to in the array
      from_array: array of available values of which the nearest has to be 
      found
   Returns:
      The nearest value found in the array
      The difference between the requested and available closest value 
      in the array
   """
   from numpy import ones,argmin
   deltas = ones(len(from_array))*1000
   for v,i in zip(from_array,range(len(from_array))):
      deltas[i] = abs(to_point - v)

   return from_array[argmin(deltas)]

def regrid_df(df,variables = [],resolution=[0]):
    import numpy as np
    from scipy.interpolate import griddata
    import pandas as pd

    string_columns = df.select_dtypes(include=['object'])

    if not len(resolution)==2:
        grid_y, grid_x = np.mgrid[
                    df['y'].min() : df['y'].max() : resolution[0],
                    df['x'].min() : df['x'].max() : resolution[0],
                    ]
    else:
        grid_y, grid_x = np.mgrid[
                    df['y'].min() : df['y'].max() : resolution[1]*1j,
                    df['x'].min() : df['x'].max() : resolution[0]*1j,
                    ]
        
    df_interpolated = pd.DataFrame({
            'x' : grid_x.ravel(),
            'y' : grid_y.ravel(),
            })

    if not len(variables):
        variables = [f for f in df.columns if not f == 'x' and not f == 'y']
    for v in variables:
        grid_var = griddata(
                ( df['x'].values , df['y'].values) , 
                df[v].values, 
                (grid_x,grid_y),
                method='linear'
                )

        df_interpolated[v] = grid_var.ravel()
        df_interpolated = df_interpolated.fillna(0)
        
    # Re-center the array to the TE location at (0,0) ##########################
    df_interpolated.y = df_interpolated.y - \
            find_nearest(0,df_interpolated.y.values)
    df_interpolated.x = df_interpolated.x - \
            find_nearest(0,df_interpolated.x.values)
    ############################################################################

    for col in string_columns.columns:
        df_interpolated[col] = df[col].unique()[0]

    return df_interpolated

def rotate_df(df, degrees = 0):
    from math import radians
    from numpy import sin, cos
    from copy import copy

    angle = radians(degrees)

    x = copy(df['x'].values)
    y = copy(df['y'].values)
    u = copy(df['u'].values)
    v = copy(df['v'].values)

    # Coordinates ##############################################################
    df['x'] =  x*cos(angle) + y*sin(angle) 
    df['y'] = -x*sin(angle) + y*cos(angle) 
    ############################################################################

    # Velocity components ######################################################
    df['u'] =  u*cos(angle) + v*sin(angle)
    df['v'] = -u*sin(angle) + v*cos(angle)
    ############################################################################

    return df


def correct_flow_plane_df( df, 
                          rotation_angle                   = 0,
                          height_correction                = 0,
                          streamwise_correction            = 0,
                         ):
    # Do all the plane correcions ########################################
    if rotation_angle:
        df = rotate_df( 
            df                               = df,
            degrees                          = rotation_angle,
        )
    ######################################################################

    # Correct height and ignore below minimum trustable y ################
    df.y = df.y + height_correction
    ######################################################################

    # Correct streanwise translation #####################################
    df.x = df.x - streamwise_correction
    ######################################################################

    return df

def rename_df_columns_from_DaVis_to_standard(df):
    DaVis_naming_dict= {
          "x"  : "x",
          "y"  : "y",
          "z"  : "z",
          "Vx" : "u",
          "Vy" : "v",
          "Vz" : "w",
          }

    df.columns = [
        DaVis_naming_dict[col] for col in df.columns
    ]

    return df

def raw_data_to_pandas_hdf5( 
    case            = '',
    root            = '' ,
    output_root     = '',
    overwrite       = False,
    time_step_limit = 0,
    plot            = False,
    airfoil_normal  = False,
):
    from os.path import join
    import case_dict_overall_correction as case_dict

    case_constants = case_dict.return_case_df()

    if not root:
        root = join(
            '/media/carlos/6E34D2CD34D29783/2015-02_SerrationPIV',
            'TR_Data_NewProcessing'
        )

    if not output_root:
        output_root = join(
            '/media/carlos/6E34D2CD34D29783/2015-02_SerrationPIV',
            'TR_Data_Location_Calibrated'
        )

    if not case: # Just to give a default value
        case = 'Sr20R21_a0_p0_U20_z05_tr_New.hdf5'

    mask = return_mask(case)

    # Get the correction values for this specific case #########################
    angle_correction = case_constants[
        case_constants.file == case
    ].rotation.values[0]

    height_correction = case_constants[
        case_constants.file == case
    ].height_correction.values[0]

    streamwise_correction = case_constants[
        case_constants.file == case
    ].x_corr.values[0]
    ############################################################################

    print "   Running case:\n{0}".format(case)
    if '.hdf5' in case:
        read_raw_hdf5_case_and_write_pandas_hdf5(
            hdf5_file             = case,
            root                  = root,
            serration_angle       = 0, # It's already corrected in the HDF5
            angle_correction      = angle_correction,
            height_correction     = height_correction,
            streamwise_correction = streamwise_correction,
            overwrite             = overwrite,
            time_step_limit       = time_step_limit,
            plot                  = plot,
            output_file           = join(output_root,case),
            airfoil_normal        = airfoil_normal,
        )
    else:
        read_raw_tecplot_case_and_write_pandas_hdf5(
            case_folder           = case,
            root                  = root,
            serration_angle       = 0, # It's already corrected in the HDF5
            angle_correction      = angle_correction,
            height_correction     = mask[1][0] + height_correction,
            streamwise_correction = mask[1][1] + streamwise_correction,
            overwrite             = overwrite,
            time_step_limit       = time_step_limit,
            output_file           = join(output_root,case),
            airfoil_normal        = airfoil_normal,
        )

def read_raw_tecplot_case_and_write_pandas_hdf5(
    case_folder,
    root                  = 0,
    output_file           = 0,
    serration_angle       = 0,
    angle_correction      = 0,
    height_correction     = 0,
    streamwise_correction = 0,
    overwrite             = False,
    time_step_limit       = 0,
    airfoil_normal        = False,
):
    from os.path import isfile,join,splitext
    from os import listdir
    from progressbar import ProgressBar,Percentage,Bar,ETA,SimpleProgress
    from pandas import HDFStore

    # File related things ######################################################
    if not output_file:
        output_file = case_folder+".hdf5"

    if airfoil_normal:
        output_file = output_file+"_AirfoilNormal"

    if not output_file.endswith('.hdf5'):
        output_file = output_file.replace(".hdf5","")+".hdf5"

    if isfile(output_file) and not overwrite:
        print "  Exiting; file exists:\n{0}".format(output_file)
        return 0
    else:
        print "  Writing\n{0}".format(output_file)
    # ##########################################################################


    hdf = HDFStore(output_file)

    time_step_files = sorted([f for f in listdir(join(root,case_folder)) \
             if splitext(f)[1] == '.dat'])

    if time_step_limit:
        time_step_files = time_step_files[:time_step_limit]

    progress = ProgressBar(
         widgets=[
             Bar(),' ',
             Percentage(),' ',
             ETA(), ' (file ',
             SimpleProgress(),')'], 
         maxval=len(time_step_files)
         ).start()

    cnt = 0
    for f,t in zip(time_step_files,range(len(time_step_files))):

       df_t = read_tecplot_file_and_correct_for_location_rotation(
           tecplot_file          = join(root,case_folder,f),
           serration_angle       = serration_angle,
           angle_correction      = angle_correction,
           height_correction     = height_correction,
           streamwise_correction = streamwise_correction,
           time_step             = t,
           airfoil_normal        = airfoil_normal,
       )

       df_t = get_vorticity(df_t)

       if cnt == 0:
           df = df_t.copy()
       else:
           df = df.append( df_t, drop_index = True)
           #df = df.drop_duplicates()

           try:
               x_cnt = df.x.value_counts().max() 
           except AttributeError:
               print df
               raise
           if not x_cnt.max() == x_cnt.min():
               print "  There's something wrong, counted {0} instances of x"\
                       .format(x_cnt.max())
               return 0

       if t == 30:
           hdf.put(case_folder,
                   df.convert_objects(), 
                   format='table', data_columns=True
                  )
       elif cnt == 30 and not t == cnt:
           hdf.append(case_folder,
                      df.convert_objects(), 
                      format='table', data_columns=True
                     )
           cnt = 0

       cnt += 1

       progress.update(t)

    progress.finish()
    hdf.close()

    return 1

def return_mask(case):
    import Masks as masks

    # Device, phi, alpha
    device,phi,alpha,U,loc,reprocessed = \
            get_case_details_from_filename(case)
    alpha = float(alpha)
    phi   = float(phi)

    # Mask
    mask_name = "{0}_phi{1:d}_alpha{2:d}_U{3}_loc{4}.dat"\
            .format(device,int(phi),int(alpha),U,loc)
    mask = masks.Masks[mask_name]

    return mask

def get_serration_angle(case):
    from math import atan,degrees

    mask = return_mask(case)

    angle = atan( 
        (mask[2][0] - mask[1][0]) / (mask[2][1] - mask[1][1])
    )
    return degrees(angle)

def read_tecplot_file_and_correct_for_location_rotation(
    tecplot_file,
    serration_angle       = 0,
    angle_correction      = 0,
    height_correction     = 0,
    streamwise_correction = 0,
    time_step             = 0,
    airfoil_normal        = False,
):

    """Reads in a tecplot file, given, and returns a pandas data frame

    Important!
    This data frame that is returned is already CORRECTED and turned to
    the standard coordinate system

    Input: address to tecplot_file

    Output: a data frame containing the 0.25 mm resolved structured grid data

    """
    import pandas as pd 
    from re import findall
    from copy import copy

    # Get available variables
    f = open(tecplot_file,'ro')

    # Read the header and get variable and title information from it ###########
    var_flag = False
    dev_flag = False
    for line in f:
        string = findall("^VARIABLES[ _A-Za-z0-9,\"=]+",line)
        if string:
            variables = [v.replace(' ','_').replace("\"","") \
                         for v in string[0].replace("VARIABLES = ",'')\
                         .split(", ")]
            variables = [v for v in variables if len(v)]
            var_flag = True
        string = findall("^TITLE = [ -_A-Za-z0-9,\"=]+",line)
        if string:
            dev_flag = True
        if var_flag and dev_flag:
            break
    f.close()
    ############################################################################

    # Put the data into a data frame ###########################################
    df = pd.read_table(
            tecplot_file,
            skiprows  = 4,
            names     = variables,
            sep       = '[ \t]+',
            index_col = False
            )
    ############################################################################
    df = df.drop('z',1)
    df = rename_df_columns_from_DaVis_to_standard(df)

    # Put the coordinate system in the standar direction, not vertical #########
    x = copy(df.x.values)
    y = copy(df.y.values)
    u = copy(df.u.values)
    v = copy(df.v.values)
    df.x =  y
    df.y = -x
    df.u =  v
    df.v = -u
    ############################################################################

    #len_x, len_y = len(df.x.unique()), len(df.y.unique())

    rotation_angle = serration_angle + angle_correction

    if airfoil_normal:
        rotation_angle = rotation_angle - 11.4

    df = correct_flow_plane_df(
        df,
        rotation_angle                   = rotation_angle,
        height_correction                = height_correction,
        streamwise_correction            = streamwise_correction,
    )

    df = regrid_df(
        df, 
        #resolution = ( len_x , len_y )
        resolution = [0.5]
    )

    df[ 'time_step' ] = time_step

    return df

#def get_time_resolved_wall_normal_line(
#    hdf5_file, case, 
#    x_locs                = [0],
#    output_hdf            = 'data.hdf5',
#    plot                  = False,
#    overwrite             = False,
#    height_correction     = 0,
#    streamwise_correction = 0,
#    rotation_correction   = 0,
#):
#   """ Gets a wall normal line, interpolated data and creates a 
#   time-resolved data frame with its information
#
#   Input:
#       hdf5_file: HDF5 file to read from
#       case: case name to read from
#       x: the 2h normalized location to read, where x=0 is the airfoil TE
#       output_hdf: how to name the pickle of the resulting data frame
#   """
#
#   #######################################################
#   #######################################################
#   # IMPORTANT
#   #
#   # The coordinates coming from the HDF5 file are the
#   # vertical freestream coordinates of DaVis.
#   #
#   # The coordinates used for the local variables are
#   # already put to the left-to-right freestream 
#   # coordinates
#   #
#   # Further, all coordinates are normalized to the
#   # tooth length, 2h
#   #
#   #######################################################
#   #######################################################
#
#   from progressbar import ProgressBar,Percentage,Bar,ETA,SimpleProgress
#   import h5py
#   import numpy as np
#   import pandas as pd
#   from scipy.interpolate import griddata
#   from copy import copy
#   from Functions import find_nearest
#   from os.path import isfile
#
#   if isfile(output_hdf) and not overwrite:
#       print "   Pickle already exists; skipping, or specify "\
#               +"overwrite next time"
#       return 0
#
#   print "   Extracting locations {0} \n".format(x_locs)\
#           +"   for case {0} from its HDF5 database\n".format(case)\
#           +"   into HDF {0}".format(output_hdf)
#
#   tooth_length = 40.
#
#   aq_freq = 10000.
#
#
#   wall_normal_lines_DF = pd.DataFrame(
#       columns = [
#           'x', 
#           'y',
#           'vx',
#           'vy',
#           'vz'
#       ]
#   )
#
#   progress = ProgressBar(
#        widgets=[
#            Bar(),' ',
#            Percentage(),' ',
#            ETA(), ' (time step ',
#            SimpleProgress(),')'], 
#        maxval=len(available_times) 
#        ).start()
#
#   t_x_cnt = 0
#   cnt = 50
#
#   df_columns = [
#       'x' ,
#       'y' ,
#       'vx',
#       'vy',
#       'vz'    ,
#   ]
#
#   hdf = pd.HDFStore(output_hdf)
#   df_xloc_wall_normal = pd.DataFrame( columns = df_columns )
#
#   for ti in available_times:
#       df['vx'] =  \
#               np.array(h5["{0}/{1}/{2}".format(case,ti,'Vy')].value)
#       df['vy'] = \
#               -np.array(h5["{0}/{1}/{2}".format(case,ti,'Vx')].value)
#       df['vz'] =  \
#               np.array(h5["{0}/{1}/{2}".format(case,ti,'Vz')].value)
#
#       df['vx_rot'],df['vy_rot'] = rotate(
#           df.vx, df.vy, flow_angle
#       )
#
#       df['x_rot'],df['y_rot'] = rotate(
#           df.x, df.y, angle
#       )
#
#       grid_x, grid_y = np.mgrid[
#                   df['x_rot'].min():df['x_rot'].max():150j,
#                   df['y_rot'].min():df['y_rot'].max():75j
#                   ]
#       grid_vx = griddata(
#               (df['x_rot'].values,df['y_rot'].values), 
#               df['vx_rot'].values, 
#               (grid_x,grid_y),
#               method='cubic'
#               )
#       grid_vy = griddata(
#               (df['x_rot'].values,df['y_rot'].values), 
#               df['vy_rot'].values, 
#               (grid_x,grid_y),
#               method='cubic'
#               )
#       grid_vz = griddata(
#               (df['x_rot'].values,df['y_rot'].values), 
#               df['vz'].values, 
#               (grid_x,grid_y),
#               method='cubic'
#               )
#
#       df_interpolated = pd.DataFrame({
#               'x' : grid_x.ravel(),
#               'y' : grid_y.ravel(),
#               'vx': grid_vx.ravel(),
#               'vy': grid_vy.ravel(),
#               'vz': grid_vz.ravel(),
#               })
#
#       # Re-center the array to the TE location at (0,0)
#       df_interpolated.y = df_interpolated.y - \
#               find_nearest(0,df_interpolated.y.values)[0]
#       df_interpolated.x = df_interpolated.x - \
#               find_nearest(0,df_interpolated.x.values)[0]
#
#       # Rotate x_locs
#       x_locs_rotated = []
#       for x_loc in x_locs:
#           x_loc_rotated,y0 = rotate(x_loc,0,angle)
#           x_locs_rotated.append(x_loc_rotated)
#           
#       if plot and ti == 0:
#           quick_plot_time_step(df       = df_interpolated,
#                                mask     = mask_rot,
#                                variable = 'vx',
#                                x_locs   = x_locs_rotated
#                               )
#
#       for x_loc_rotated in x_locs_rotated:
#
#           nearest_available_x,tmp = find_nearest(
#               x_loc_rotated,
#               df_interpolated[df_interpolated.y==0].x.values
#           )
#
#           df_xloc_wall_normal_step = df_interpolated[
#               (df_interpolated.x == nearest_available_x) &\
#               (df_interpolated.y > 0) & \
#               (df_interpolated.y < 0.95*df_interpolated.y.max()) 
#           ]
#
#           df_xloc_wall_normal_step['ti'] = ti
#           df_xloc_wall_normal_step['t_real'] = ti/aq_freq
#
#           if df_xloc_wall_normal_step.empty:
#               print "Error, available rotated x coordinate not in mesh"
#               h5.close()
#               progress.finish()
#               print x_loc_rotated
#               break
#
#
#           if cnt == 50:
#               if t_x_cnt == cnt:
#                   hdf.put(case,
#                           df_xloc_wall_normal.convert_objects(), 
#                           format='table', data_columns=True
#                          )
#               else:
#                   hdf.append(case,
#                              df_xloc_wall_normal.convert_objects(), 
#                              format='table', data_columns=True
#                             )
#
#               df_xloc_wall_normal = pd.DataFrame( 
#                   columns = df_columns 
#               )
#               cnt = 0
#           else:
#               df_xloc_wall_normal = df_xloc_wall_normal.append(
#                   df_xloc_wall_normal_step, ignore_index=True
#               )
#               cnt += 1
#
#       t_x_cnt += 1
#       progress.update(t_x_cnt)
#
#   hdf.close()
#   h5.close()
#
#   progress.finish()
#   return wall_normal_lines_DF
#

def run_test(df):
    import numpy as np
    import matplotlib.pyplot as plt

    y_locs = df.y.unique()

    u_avg = np.zeros(len(y_locs))
    for y,i in zip(y_locs,range(len(y_locs))):
        u_avg[i] = df[df.y == y].u.mean()

    plt.plot( u_avg, y_locs )
    plt.show()

def plot_first_frame_of_hdf(hdf_file):
    import pandas as pd
    from os.path import split, splitext
    
    df = pd.read_hdf( hdf_file,
                     splitext(split(hdf_file)[1])[0]\
                     .replace("_AirfoilNormal",""),
                     where = [ 'time_step = 0' ],
                     columns = ['x']
                    )
    print df.x.values_count()
    #show_surface_from_df(df)
    print "OK: {0}".format(hdf_file)

def average_time_series(df):
    import pandas as pd
    from progressbar import ProgressBar,Percentage,Bar,ETA,SimpleProgress

    averaged_df = pd.DataFrame(columns = ['x','y','u','v','w',
                                'u_rms','v_rms','w_rms'])

    grouped_df = df.groupby(['x','y'])

    progress = ProgressBar(
         widgets=[
             Bar(),' ',
             Percentage(),' ',
             ETA(), ' (file ',
             SimpleProgress(),')'], 
         maxval=len(grouped_df)
         ).start()

    cnt = 0
    for (x, y), data in grouped_df:
        result = {
            'x':x,
            'y':y,
        }
        for v in ['u','v','w']:
            result[v]        = data[v].mean()
            result[v+'_rms'] = data[v].std()

        averaged_df = averaged_df.append( pd.DataFrame( result , index = [0]), 
                                         ignore_index = True)
        cnt += 1
        progress.update(cnt)

    progress.finish()
    return averaged_df

def get_point_time_series_from_pandas_hdf( hdf, x,y , overwrite = False):
    import pandas as pd
    import numpy as np
    from os.path import isfile, join, split

    time_series_pickle_name = "{0}_x{1}_y{2}.p".format(
        split(hdf)[1].replace('.hdf',''), x, y
    )

    if isfile( join( "ReservedData",time_series_pickle_name ) ) \
       and not overwrite:
        time_series = pd.read_pickle( 
            join( "ReservedData",time_series_pickle_name ) 
        )
        return time_series

    # First we need to find out what coordinates are available #################
    coords = pd.read_hdf(
        hdf,
        split(hdf)[1].replace('.hdf5','').replace('_AirfoilNormal',''),
        where = ['time_step = 0'],
        columns = [ 'x', 'y' ],
    )
    # ##########################################################################

    # Now we need to find the coordinate pair that closest gets to the #########
    # requested one ############################################################
    diff_x = np.abs( x - coords.x )
    diff_y = np.abs( y - coords.y )
    selected_x = coords.x.values[np.argmin(diff_x + diff_y)]
    selected_y = coords.y.values[np.argmin(diff_x + diff_y)]

    if abs(selected_x - x) > 1.:
        print " Found a too large distance from the requested x"
        print " Requested {0}".format(x)
        print " Closest available {0}".format(selected_x)
        return 0
    if abs(selected_y - y) > 1.:
        print " Found a too large distance from the requested y"
        print " Requested {0}".format(y)
        print " Closest available {0}".format(selected_y)
        return 0

    print " Found the following closest points:"
    print " Asked x = {0}".format(x)
    print " Found x = {0}".format(selected_x)
    print " Asked y = {0}".format(y)
    print " Found y = {0}".format(selected_y)
    # ##########################################################################

    # Now get the time series ##################################################
    time_series = pd.read_hdf(
        hdf,
        split(hdf)[1].replace('.hdf5','').replace('_AirfoilNormal',''),
        where = [ 
            'x < {0}'.format( selected_x + 0.01),
            'x > {0}'.format( selected_x - 0.01),
            'y < {0}'.format( selected_y + 0.01),
            'y > {0}'.format( selected_y - 0.01),
                ],
        columns = ['u','v','w','time_step','x','y']
    )
    # ##########################################################################

    time_series = time_series.sort_values( by = 'time_step' )
    time_series.to_pickle(
        join( 'ReservedData', '{0}'.format(time_series_pickle_name) )
    )
    return time_series, coords.x, coords.y

def wall_normal_data_to_reserved_pickles_from_pandas_hdf( hdf, x ):
    from numpy import array

    y_locs = array([0.1, 0.3, 0.5, 0.6, 1.0, 1.5]) * 9.

    for y in y_locs:
        get_point_time_series_from_pandas_hdf( hdf, x, y, overwrite = False)




