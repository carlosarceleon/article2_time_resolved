def build_data_hdf5(root,case_folders,target,overwrite=False):
   """ Builds an HDF5 file containing the given cases' data time series

   Input:
        root: root location of the case folders
        case_folders: name of the case folders to process as a list
        target: the HDF5 to write
        overwrite: if an HDF5 already exists, overwrite it (bool)
   """
   
   from progressbar import ProgressBar,Percentage,Bar,ETA,SimpleProgress
   from Functions import read_tecplot_file, get_case_details_from_filename
   import Masks as masks
   import h5py
   import os
   from math import atan
   from numpy import deg2rad

      
   acquisition_frequency = 1./5000.

   # Find number of total files I need to process ##############################
   print "   Going to process files from the folders:"
   n_files = 0
   for cf in [case_folders]:
      n_files += len([f for f in os.listdir(os.path.join(root,cf)) \
                      if f.endswith('.dat')])
      print "      {0}".format(cf)
   #############################################################################

   # Check if the file already exists, otherwise start writing #################
   if os.path.isfile(target):
      if os.path.getsize(target) < 10000 or overwrite:
         os.remove(target)
      else:
         print "    File exists, not overwriting\n"
         return 1
   print "  Saving to {0}".format(target)

   try:
       h5 = h5py.File(target+'.hdf5','w')
   except:
       return 0
   #############################################################################

   progress = ProgressBar(
        widgets=[
            Bar(),' ',
            Percentage(),' ',
            ETA(), ' (file ',
            SimpleProgress(),')'], 
        maxval=n_files
        ).start()

   # Run through all folders ###################################################
   cnt_files = 0
   for cf in [case_folders]:

      # Run through all time step datafiles that were found in the folder ######
      files = [f for f in os.listdir(os.path.join(root,cf)) \
               if os.path.splitext(f)[1] == '.dat']

      for f,t in zip(files,range(len(files))):

         # If it's the first time step, initialize the hdf5 group ##############
         df = read_tecplot_file(os.path.join(root,cf,f))
         if f == files[0]:
            grp = h5.create_group(cf)

            # Coordinate points (number of)
            planar_data=False

            grp.attrs['nx'] = df.x.size
            grp.attrs['ny'] = df.y.size

            try:
                grp.attrs['nz'] = df.z.size
            except AttributeError:
                planar_data = True
            
            # Device, phi, alpha
            device,phi,alpha,U,loc,reprocessed = \
                    get_case_details_from_filename(cf)
            alpha = float(alpha)
            phi = float(phi)

            # Mask
            mask_name = "{0}_phi{1:d}_alpha{2:d}_U{3}_loc{4}.dat"\
                    .format(device,int(phi),int(alpha),U,loc)
            mask = masks.Masks[mask_name]

            # Rotation angle so that true Vy is vertical (and streamwise)
            if alpha: sign = alpha/abs(alpha)
            else: sign = 1
            if alpha == -6:
                alpha = -12
            angle = atan( 
                (mask[2][0] - mask[1][0]) / (mask[2][1] - mask[1][1])
            )
            grp.attrs['mask_name']    = mask_name
            grp.attrs['device']       = device
            grp.attrs['phi']          = phi
            grp.attrs['alpha']        = alpha
            grp.attrs['U_inf']        = U
            grp.attrs['loc']          = loc
            grp.create_dataset('mask', data=mask)
            grp.attrs['angle']        = angle
            grp.attrs['flow_angle']   = angle + sign \
                    * deg2rad(abs(phi)+abs(alpha))
            # Coordinate points 
            grp.create_dataset('x', 
                               data = df.x.values-masks.Masks[mask_name][1][0],
                               dtype ='float')
            grp.create_dataset('y', 
                               data = df.y.values-masks.Masks[mask_name][1][1],
                               dtype ='float')
      
         # Create a new group to store the datasets for this time
         grp = h5.create_group("{0}/{1}".format(cf,t))
         grp.attrs['time'] = t*acquisition_frequency
         grp.create_dataset('Vx', data= df['Vx'].values,dtype='float')
         grp.create_dataset('Vy', data= df['Vy'].values,dtype='float')
         if not planar_data:
             grp.create_dataset('Vz', data= df['Vz'].values,dtype='float')

            
         cnt_files+=1
         progress.update(cnt_files)

   progress.finish()

   h5.close()

def read_hdf5_time_series(hdf5_file,case,loc=(-0.1,1)):
   """ Reads an HDF5 file and returns the velocity components time series

   Input:
        hdf5_file: which HDF5 file to read
        case: which case to extract the data from
        (x,y): normalized coordinate location where to extract the data 
            from. (0,0) is the TE. This is in DaVis coordinates,
            meaning that -x is the airfoil y coordinate

   Output:
        df: DataFrame
        Returns zero if the coordinate is beyond the frame of view
   """


   from progressbar import ProgressBar,Percentage,Bar,ETA,SimpleProgress
   import h5py
   import numpy as np
   import pandas as pd
   import re
   from Functions import find_nearest

   x = loc[0]; y = loc[1]

   tooth_length = 40.

   h5 = h5py.File(hdf5_file,'r')
   # Read all the available time steps
   tn = sorted(map(int,re.findall("[0-9]+"," ".join([n for n in h5[case]]))))
   # Prepare the velocity arrays
   vx = np.zeros(len(tn))
   vy = np.zeros(len(tn))
   vz = np.zeros(len(tn))
   time  = np.zeros(len(tn))

   # Read the available coordinates
   h5_x = h5["{0}/x".format(case)].value/tooth_length
   h5_y = h5["{0}/y".format(case)].value/tooth_length


   # This is just 'angle' so that the flow is analized
   # as going parallel to the serration surface, and not
   # in the direction of the standard coordinate system
   angle = -h5["{0}".format(case)].attrs['angle'] 

   x_rot, y_rot = rotate(x,y,angle)

   xi,xn,xd = find_nearest(x_rot,np.unique(h5_x))
   yi,yn,yd = find_nearest(y_rot,np.unique(h5_y))

   if x_rot<min(h5_x) or x_rot>max(h5_x)\
      or y_rot<min(h5_y) or y_rot>max(h5_y):
       return None

   ind = (xn==h5_x)*(yn==h5_y) # TODO: IS THIS RIGHT?!

   print "   Looking for coordinate points close to x = {0:.2f}"\
           +" and y = {1:.2f}".format(x_rot,y_rot)

   print "      Found x = {0:.2f} and y = {1:.2f}".format(
       float(h5_x[ind]),float(h5_y[ind])
   )

   progress = ProgressBar(
        widgets=[
            Bar(),' ',
            Percentage(),' ',
            ETA(), ' (timestep ',
            SimpleProgress(),')'], 
        maxval=len(tn)
        ).start()
   for t,i in zip(tn,range(len(tn))):
      vx[i] = h5["{0}/{1}/{2}".format(case,t,'Vx')][ind]
      vy[i] = h5["{0}/{1}/{2}".format(case,t,'Vy')][ind]
      vz[i] = h5["{0}/{1}/{2}".format(case,t,'Vz')][ind]
      time[i]  = h5["{0}/{1}".format(case,t)].attrs['time']
      progress.update(i)
   progress.finish()
   h5.close()

   print "   Rotating to angle {0:.2f}".format(np.rad2deg(angle))
   vx_rot, vy_rot = rotate(vx,vy,angle)
      
   return pd.DataFrame({"vx":vx_rot,'vy':vy_rot,"vz":vz,"t":time})

def rotate(x,y,angle):
    from math import sin,cos
    from numpy import array

    x = array(x)
    y = array(y)

    x_rot =  x*cos(angle) + y*sin(angle)
    y_rot = -x*sin(angle) + y*cos(angle)

    return x_rot,y_rot

def return_rotated_hdf5_in_df(
    hdf5_file,
    case,
    output_name,
    overwrite=False):
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
   from os.path import isfile

   if isfile(output_name) and not overwrite:
       print "   Pickle already exists; skipping, or specify "\
               +"overwrite next time"
       return 0

   tooth_length = 40.


   h5 = h5py.File(hdf5_file,'r')

   # Read the available coordinates
   
   x =  np.array(h5["{0}/y".format(case)].value)/tooth_length
   y = -np.array(h5["{0}/x".format(case)].value)/tooth_length

   # Rotate requested coordinate points
   angle = -h5["{0}".format(case)].attrs['angle']

   flow_angle = angle

   x_rot, y_rot = rotate(
       x, y, angle
   )
   
   available_times = sorted([int(f[0]) for f in \
                             h5['{0}'.format(case)].iteritems()\
                             if not 'mask' in f and not 'x' in f and not 'y'\
                             in f])

   progress = ProgressBar(
        widgets=[
            Bar(),' ',
            Percentage(),' ',
            ETA(), ' (time step ',
            SimpleProgress(),')'], 
        maxval=len(available_times)
        ).start()

   t_x_cnt = 0
   df_columns = [
       'x_rot' ,
       'y_rot' ,
       'vx_rot',
       'vy_rot',
       'vz'    ,
       'ti'    ,
   ]
   df_ti = pd.DataFrame( columns = df_columns )
   cnt = 0
   hdf = pd.HDFStore(output_name)
   for ti in available_times[:500]:

       df_ti_step = pd.DataFrame( data = {
           'x_rot' : x_rot,
           'y_rot' : y_rot,
       })

       vx               = np.array(h5["{0}/{1}/{2}".format(case,ti,'Vy')]\
                                   .value)
       vy               =-np.array(h5["{0}/{1}/{2}".format(case,ti,'Vx')]\
                                   .value)
       df_ti_step['vz'] = np.array(h5["{0}/{1}/{2}".format(case,ti,'Vz')]\
                                   .value)

       df_ti_step['vx_rot'],df_ti_step['vy_rot'] = rotate(
           vx, vy, flow_angle
       )

       df_ti_step['ti']     = ti

       t_x_cnt += 1
       cnt     += 1

       progress.update(t_x_cnt)

       if cnt == 50:
           if t_x_cnt == cnt:
               hdf.put(case    , df_ti.convert_objects(), 
                       format='table', data_columns=True)
           else:
               hdf.append(case , df_ti.convert_objects(), 
                          format='table', data_columns=True)

           df_ti = pd.DataFrame( columns = df_columns)
           cnt = 0
       else:
           df_ti = df_ti.append(df_ti_step,ignore_index=True)

   hdf.close()
   h5.close()

   progress.finish()

def rotate_df(df, degrees = 0):
    from math import radians
    from numpy import sin, cos

    angle = radians(degrees)

    x = df['x'] 
    y = df['y'] 
    u = df['u']
    v = df['v']

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
                          rotation_angle        = 0,
                          height_correction     = 0,
                          streamwise_correction = 0,
                         ):
    # Do all the plane correcions ########################################
    if rotation_angle:
        df = rotate_df( df, rotation_angle )
    ######################################################################

    # Correct height and ignore below minimum trustable y ################
    df.y = df.y+height_correction
    ######################################################################

    # Correct streanwise translation #####################################
    df.x = df.x-streamwise_correction
    ######################################################################

    return df

def get_time_resolved_wall_normal_line(
    hdf5_file, case, 
    x_locs                = [0],
    output_hdf            = 'data.hdf5',
    plot                  = False,
    overwrite             = False,
    height_correction     = 0,
    streamwise_correction = 0,
    rotation_correction   = 0,
):
   """ Gets a wall normal line, interpolated data and creates a 
   time-resolved data frame with its information

   Input:
       hdf5_file: HDF5 file to read from
       case: case name to read from
       x: the 2h normalized location to read, where x=0 is the airfoil TE
       output_hdf: how to name the pickle of the resulting data frame
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
   # Further, all coordinates are normalized to the
   # tooth length, 2h
   #
   #######################################################
   #######################################################

   from progressbar import ProgressBar,Percentage,Bar,ETA,SimpleProgress
   import h5py
   import numpy as np
   import pandas as pd
   from scipy.interpolate import griddata
   from copy import copy
   from Functions import find_nearest
   from os.path import isfile

   if isfile(output_hdf) and not overwrite:
       print "   Pickle already exists; skipping, or specify "\
               +"overwrite next time"
       return 0

   print "   Extracting locations {0} \n".format(x_locs)\
           +"   for case {0} from its HDF5 database\n".format(case)\
           +"   into HDF {0}".format(output_hdf)

   tooth_length = 40.

   aq_freq = 10000.

   h5 = h5py.File(hdf5_file,'r')

   # Read the available coordinates
   df = pd.DataFrame({
       "x": np.array(h5["{0}/y".format(case)].value)/tooth_length,
       "y":-np.array(h5["{0}/x".format(case)].value)/tooth_length
       })

   # Rotate requested coordinate points
   angle = -h5["{0}".format(case)].attrs['angle']

   flow_angle = angle

   mask = h5["{0}/mask".format(case)].value/40.
   mask = mask - mask[1]

   mask_rot = copy(mask)

   for mi in range(len(mask)):
       mask_rot[mi][0], mask_rot[mi][1] = rotate(
           mask[mi][0],mask[mi][1],angle
       )

   df['x_rot'],df['y_rot'] = rotate(
       df.x, df.y, angle
   )
   
   available_times = sorted(
       [int(f[0]) for f in \
        h5['{0}'.format(case)].iteritems()\
        if not 'mask' in f and not 'x' in f and not 'y'\
        in f]
   )

   wall_normal_lines_DF = pd.DataFrame(
       columns = [
           'x', 
           'y',
           'vx',
           'vy',
           'vz'
       ]
   )

   progress = ProgressBar(
        widgets=[
            Bar(),' ',
            Percentage(),' ',
            ETA(), ' (time step ',
            SimpleProgress(),')'], 
        maxval=len(available_times) 
        ).start()

   t_x_cnt = 0
   cnt = 50

   df_columns = [
       'x' ,
       'y' ,
       'vx',
       'vy',
       'vz'    ,
   ]

   hdf = pd.HDFStore(output_hdf)
   df_xloc_wall_normal = pd.DataFrame( columns = df_columns )

   for ti in available_times:
       df['vx'] =  \
               np.array(h5["{0}/{1}/{2}".format(case,ti,'Vy')].value)
       df['vy'] = \
               -np.array(h5["{0}/{1}/{2}".format(case,ti,'Vx')].value)
       df['vz'] =  \
               np.array(h5["{0}/{1}/{2}".format(case,ti,'Vz')].value)

       df['vx_rot'],df['vy_rot'] = rotate(
           df.vx, df.vy, flow_angle
       )

       df['x_rot'],df['y_rot'] = rotate(
           df.x, df.y, angle
       )

       grid_x, grid_y = np.mgrid[
                   df['x_rot'].min():df['x_rot'].max():150j,
                   df['y_rot'].min():df['y_rot'].max():75j
                   ]
       grid_vx = griddata(
               (df['x_rot'].values,df['y_rot'].values), 
               df['vx_rot'].values, 
               (grid_x,grid_y),
               method='cubic'
               )
       grid_vy = griddata(
               (df['x_rot'].values,df['y_rot'].values), 
               df['vy_rot'].values, 
               (grid_x,grid_y),
               method='cubic'
               )
       grid_vz = griddata(
               (df['x_rot'].values,df['y_rot'].values), 
               df['vz'].values, 
               (grid_x,grid_y),
               method='cubic'
               )

       df_interpolated = pd.DataFrame({
               'x' : grid_x.ravel(),
               'y' : grid_y.ravel(),
               'vx': grid_vx.ravel(),
               'vy': grid_vy.ravel(),
               'vz': grid_vz.ravel(),
               })

       # Re-center the array to the TE location at (0,0)
       df_interpolated.y = df_interpolated.y - \
               find_nearest(0,df_interpolated.y.values)[0]
       df_interpolated.x = df_interpolated.x - \
               find_nearest(0,df_interpolated.x.values)[0]

       # Rotate x_locs
       x_locs_rotated = []
       for x_loc in x_locs:
           x_loc_rotated,y0 = rotate(x_loc,0,angle)
           x_locs_rotated.append(x_loc_rotated)
           
       if plot and ti == 0:
           quick_plot_time_step(df       = df_interpolated,
                                mask     = mask_rot,
                                variable = 'vx',
                                x_locs   = x_locs_rotated
                               )

       for x_loc_rotated in x_locs_rotated:

           nearest_available_x,tmp = find_nearest(
               x_loc_rotated,
               df_interpolated[df_interpolated.y==0].x.values
           )

           df_xloc_wall_normal_step = df_interpolated[
               (df_interpolated.x == nearest_available_x) &\
               (df_interpolated.y > 0) & \
               (df_interpolated.y < 0.95*df_interpolated.y.max()) 
           ]

           df_xloc_wall_normal_step['ti'] = ti
           df_xloc_wall_normal_step['t_real'] = ti/aq_freq

           if df_xloc_wall_normal_step.empty:
               print "Error, available rotated x coordinate not in mesh"
               h5.close()
               progress.finish()
               print x_loc_rotated
               break


           if cnt == 50:
               if t_x_cnt == cnt:
                   hdf.put(case,
                           df_xloc_wall_normal.convert_objects(), 
                           format='table', data_columns=True
                          )
               else:
                   hdf.append(case,
                              df_xloc_wall_normal.convert_objects(), 
                              format='table', data_columns=True
                             )

               df_xloc_wall_normal = pd.DataFrame( 
                   columns = df_columns 
               )
               cnt = 0
           else:
               df_xloc_wall_normal = df_xloc_wall_normal.append(
                   df_xloc_wall_normal_step, ignore_index=True
               )
               cnt += 1

       t_x_cnt += 1
       progress.update(t_x_cnt)

   hdf.close()
   h5.close()

   progress.finish()
   return wall_normal_lines_DF

def plot_surface(hdf5_file,case,time=[0],variable="Vx",
                 plt_name='test.png',p=(0,0)):
   """ Plots the velocity values for the given component at the given time

   Input:
       hdf5_file: HDF5 file to read from
       case: case name to read from
       time: the time step to read from
       variable: the variable to read (["Vx"],"Vy","Vz")
       plt_name: the output name of the plot (default is test.png)
   """
   import h5py
   import numpy as np
   import pandas as pd
   from Functions import find_nearest
   from matplotlib import pyplot as plt
   from scipy.interpolate import griddata

   tooth_length = 40.

   h5 = h5py.File(hdf5_file,'r')

   # Read the available coordinates
   df = pd.DataFrame({
       "x": np.array(h5["{0}/y".format(case)].value)/tooth_length,
       "y":-np.array(h5["{0}/x".format(case)].value)/tooth_length
       })
   mask = h5["{0}/mask".format(case)].value/40.
   mask = mask - mask[1]

   # Rotate requested coordinate points
   angle = -h5["{0}".format(case)].attrs['angle']
   #angle = -np.deg2rad(15)
   #flow_angle = h5["{0}".format(case)].attrs['flow_angle']
   flow_angle = angle

   for mi in range(len(mask)):
       mask[mi][0], mask[mi][1] = rotate(
           mask[mi][0],mask[mi][1],angle
       )

   for ti in time:
       df['vx'] = np.array(h5["{0}/{1}/{2}".format(case,ti,'Vy')].value)
       df['vy'] = -np.array(h5["{0}/{1}/{2}".format(case,ti,'Vx')].value)
       df['vz'] = np.array(h5["{0}/{1}/{2}".format(case,ti,'Vz')].value)


       print "   Rotating the flow to angle {0:.2f}".format(
           np.rad2deg(flow_angle)
       )
       df['vx_rot'],df['vy_rot'] = rotate(
           df.vx, df.vy, flow_angle
       )

       df['x_rot'],df['y_rot'] = rotate(
           df.x, df.y, angle
       )

       grid_x, grid_y = np.mgrid[
                   df['x_rot'].min():df['x_rot'].max():75j, 
                   df['y_rot'].min():df['y_rot'].max():150j
                   ]
       grid_vx = griddata(
               (df['x_rot'].values,df['y_rot'].values), 
               df['vx_rot'].values, 
               (grid_x,grid_y),
               method='cubic'
               )
       grid_vy = griddata(
               (df['x_rot'].values,df['y_rot'].values), 
               df['vy_rot'].values, 
               (grid_x,grid_y),
               method='cubic'
               )
       grid_vz = griddata(
               (df['x_rot'].values,df['y_rot'].values), 
               df['vz'].values, 
               (grid_x,grid_y),
               method='cubic'
               )

       df_interpolated = pd.DataFrame({
               'x' : grid_x.ravel(),
               'y' : grid_y.ravel(),
               'vx': grid_vx.ravel(),
               'vy': grid_vy.ravel(),
               'vz': grid_vz.ravel(),
               })

       # Re-center the array to the TE location at (0,0)
       df_interpolated.y = df_interpolated.y - \
               find_nearest(0,df_interpolated.y.values)[0]
       df_interpolated.x = df_interpolated.x - \
               find_nearest(0,df_interpolated.x.values)[0]

       fig = plt.figure(figsize=(12,8))
       base_cmap = plt.get_cmap("RdYlBu_r")
       plt.subplot(111,aspect=1)
       cf = plt.contourf(
               grid_x,
               grid_y,
               df_interpolated[variable].reshape(grid_x.shape),
               #levels=[-200., -160., -120.,  -80.,  -40., 40.,   80.,  120.,
            #160.,  200.],
               cmap=base_cmap
               )
       df_interpolated = df_interpolated.sort('y')
       plt.quiver( df_interpolated.x.values[::20], 
                  df_interpolated.y.values[::20], 
                  df_interpolated.vx.values[::20], 
                  df_interpolated.vy.values[::20],
                  linewidths=(1,), edgecolors=('k'), scale=400 
                 )

       clb = plt.colorbar(cf)
       clb.set_label("$\\omega_z$")

       plt.ylabel("$\\tilde x/2h$")
       plt.xlabel("$\\tilde y/2h$")

       #co = plt.contour(
       #        grid_x,
       #        grid_y,
       #        grid_z,
       #        c='k'
       #        )

       plt.fill_between(
           mask[:,1],
           min(-mask[:,0]),
           -mask[:,0],
           color='k'
       )
       #plt.scatter(0,0,s=300,color='r')

       #plt.scatter(x_rot,y_rot,s=300)
       if len(time)==1:
           plt.savefig(plt_name)
       else:
           plt.savefig("{0:05d}_{1}".format(ti,plt_name),
                       bbox_inches='tight')
       plt.cla()
   h5.close()

   return df['vx_rot'],df['vy_rot'],df['vz']

def quick_plot_time_step(df,mask=False,variable='vx',x_locs=[]):
   import numpy as np
   from matplotlib import pyplot as plt

   fig                    = plt.figure(figsize=(12,8))
   base_cmap              = plt.get_cmap("RdYlBu_r")
   plt.subplot(111,aspect = 1)
   #grid_y,grid_x          = np.meshgrid(df.y.unique(),df.x.unique())
   #cf = plt.contourf(
   #        grid_x,
   #        grid_y,
   #        df[variable].reshape(grid_x.shape),
   #        cmap=base_cmap
   #        )
   plt.quiver( df.x.values[::20], 
              df.y.values[::20], 
              df.vx.values[::20], 
              df.vy.values[::20],
              linewidths=(1,), edgecolors=('k'), scale=400 
             )

   #clb = plt.colorbar(cf)
   #clb.set_label("$\\omega_z$")

   plt.ylabel("$\\tilde x/2h$")
   plt.xlabel("$\\tilde y/2h$")

   if mask:
       plt.fill_between(
           mask[:,1],
           min(-mask[:,0]),
           -mask[:,0],
           color='k'
       )
   #plt.scatter(0,0,s=300,color='r')

   if len(x_locs):
       plt.scatter(x_locs,[0]*len(x_locs),s=300)
   plt.show()
   plt.cla()

   return 0

def read_time_series_range(hdf5_file,case,variable="Vx",ti=0,tf=100):
   """ Reads an HDF5 file and returns the velocity components time series for a given range of time

   Input:
        hdf5_file: which HDF5 file to read
        case: which case to extract the data from
        (x,y): normalized coordinate location where to extract the data from. (0,0) is the TE
        ti: initial time step
        tf: final time step

   Output:
        DataFrame
   """
   from progressbar import ProgressBar,Percentage,Bar,ETA,SimpleProgress
   import h5py
   import numpy as np
   import pandas as pd
   from math import sin,cos

   tooth_length = 40.

   h5 = h5py.File(hdf5_file,'r')

   # Read the available coordinates
   mask = h5["{0}/mask".format(case)].value/40.
   mask = mask - mask[1]

   # Rotate requested coordinate points
   flow_angle = h5["{0}".format(case)].attrs['flow_angle']

   progress = ProgressBar(
        widgets=[
            Bar(),' ',
            Percentage(),' ',
            ETA(), ' (time step ',
            SimpleProgress(),')'], 
        maxval=tf-ti
        ).start()

   for t in range(ti,tf+1):
       df_tmp = pd.DataFrame({
           "x":np.array(h5["{0}/x".format(case)].value)/tooth_length,
           "y":np.array(h5["{0}/y".format(case)].value)/tooth_length
           })
       df_tmp['vx']   = np.array(h5["{0}/{1}/{2}".format(case,t,'Vx')].value)
       df_tmp['vy']   = np.array(h5["{0}/{1}/{2}".format(case,t,'Vy')].value)
       df_tmp['vz']   = np.array(h5["{0}/{1}/{2}".format(case,t,'Vz')].value)
       df_tmp['t']    = t
       df_tmp['case'] = case

       #print "   Rotating the flow to angle {0:.2f}".format(np.rad2deg(flow_angle))
       df_tmp['vx_rot'] = df_tmp['vx']*cos(flow_angle)+df_tmp['vy']*sin(flow_angle)
       df_tmp['vy_rot'] = -df_tmp['vx']*sin(flow_angle) + df_tmp['vy']*cos(flow_angle)
       
       if t == ti:
           df = df_tmp.copy()
       else:
           df = df.append(df_tmp)

       progress.update(t-ti)
   h5.close()
   progress.finish()

   return df

