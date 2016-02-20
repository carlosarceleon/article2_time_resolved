DaVis_naming_convention_dictionary = {
      "x"                            : "x"                 ,
      "y"                            : "y"                 ,
      "AvgVx"                        : "Avg Vx"            ,
      "AvgVy"                        : "Avg Vy"            ,
      "AvgVz"                        : "Avg Vz"            ,
      "LengthofAvgV"                 : "Length of Avg V"   ,
      "StandarddeviationofVx"        : "RMS Vx"            ,
      "StandarddeviationofVy"        : "RMS Vy"            ,
      "StandarddeviationofVz"        : "RMS Vz"            ,
      "LengthofStandarddeviationofV" : "Length of RMS V"   ,
      "ReynoldstressXY"              : "Reynold stress XY" ,
      "ReynoldstressXZ"              : "Reynold stress XZ" ,
      "ReynoldstressYZ"              : "Reynold stress YZ" ,
      "ReynoldstressXX"              : "Reynold stress XX" ,
      "ReynoldstressYY"              : "Reynold stress YY" ,
      "ReynoldstressZZ"              : "Reynold stress ZZ" ,
      }

DaVis_naming_convention = [
    "x"                 ,
    "y"                 ,
    "Avg Vx"            ,
    "Avg Vy"            ,
    "Avg Vz"            ,
    "Length of Avg V"   ,
    "RMS Vx"            ,
    "RMS Vy"            ,
    "RMS Vz"            ,
    "Length of RMS V"   ,
    "Reynold stress XY" ,
    "Reynold stress XZ" ,
    "Reynold stress YZ" ,
    "Reynold stress XX" ,
    "Reynold stress YY" ,
    "Reynold stress ZZ" ,
]

def get_mask_segment_functions(mask):
    from scipy.interpolate import interp1d
    from numpy import array,zeros

    s1 = 0; s2 = 0; s3 = 0; s4 = 0

    # The iffs are in case the points are in the same location, then return
    # a zero for that line segment interpolation

    x_points = [mask[0][0],mask[1][0],mask[2][0],mask[3][0],mask[4][0]]
    y_points = [mask[0][1],mask[1][1],mask[2][1],mask[3][1],mask[4][1]]

    # Find the inflection point:
    cnt = 0
    prev = y_points[0]
    for y in y_points[1:]:
        if y < prev: 
            print "      Inflection point found at node {0}, y = {1:.2f} (starting from 0)".format(cnt,y)
            break
        else:
            cnt += 1
            prev = y

    if y_points[0] == y_points[4]:
        pass
    else:
        if y_points[0] < y_points[4]:
            y_new_points = zeros(len(y_points)+1)
            y_new_points[-1] = y_points[0]
            x_new_points = zeros(len(x_points)+1)
            x_new_points[-1] = x_points[0]
            for y,x,i in zip(y_points,x_points,range(len(y_points))):
                y_new_points[i] = y
                x_new_points[i] = x
        elif y_points[0] > y_points[4]:
            y_new_points = zeros(len(y_points)+1)
            y_new_points[0] = y_points[-1]
            x_new_points = zeros(len(x_points)+1)
            x_new_points[0] = x_points[-1]
            for y,x,i in zip(y_points,x_points,range(len(y_points))):
                y_new_points[i+1] = y
                x_new_points[i+1] = x
            cnt = cnt + 1
        y_points = y_new_points
        x_points = x_new_points


    s1 = interp1d(y_points[:cnt+1]    , x_points[:cnt+1]  )
    s2 = interp1d(y_points[:cnt-1:-1] , x_points[:cnt-1:-1])

    return s1,s2,s3,s4,cnt

def apply_mask(df,case_file):
    from copy import copy
    from pandas import DataFrame
    from numpy import nan
    import Masks

    mask = Masks.Masks[case_file]

    s1,s2,s3,s4,inflection_cnt = get_mask_segment_functions(mask)
    mask_y_points = [mask[0][1],mask[1][1],mask[2][1],mask[3][1],mask[4][1]]

    df_masked = copy(df)

    # Scan all values in y to build boolean mask
    for y in df['y'].unique():
        # Ignore all values left of segment 1
        if y<max(mask_y_points) and y>min(mask_y_points):
            for var in df.columns:
                if var == 'x' or var == 'y':
                    pass
                else:
                    try: s1(y)
                    except ValueError: print "   s1 failed; ", y,mask_y_points; raise
                    try: s2(y)
                    except ValueError: print "   s2 failed; ", y,mask_y_points; raise
                    df_masked.loc[
                            (df['x'] > float(s1(y))) & 
                            (df['x'] < float(s2(y))) & 
                            (df['y'] == y) ,
                            var] = nan
                            #df[ 
                            #        (df['x'] > float(s1(y))) & 
                            #        (df['x'] < float(s2(y))) & 
                            #        (df['y'] == y) 
                            #        ][var] * nan
    return df_masked

def find_nearest(to_point,from_array):
   """ Finds the nearest available value in a array to a given value

   Inputs:
      to_point: value to find the nearest to in the array
      from_array: array of available values of which the nearest has to be found
   Returns:
      The nearest value found in the array
      The difference between the requested and available closest value in the array
   """
   from numpy import ones,argmin
   deltas = ones(len(from_array))*1000
   for v,i in zip(from_array,range(len(from_array))):
      deltas[i] = abs(to_point - v)

   return from_array[argmin(deltas)],deltas[argmin(deltas)]

def get_variables_in_tecplot(tecplot_file):
   """ Returns the available variables in this tecplot file

   Input: tecplot file path
   Output: vector of strings
   """
   import re
   from numpy import array
   variables = []
   with open(tecplot_file, 'r') as txt:
       lines = txt.readlines()

       second_line = lines[1]

       x   = re.search('VARIABLES = "(.+?)"',second_line)
       y   = re.search('VARIABLES = "[a-z]", "(.+?)"',second_line)
       var = re.search('VARIABLES = "[a-z]", "[a-z]", "(.+?)"',second_line)
       try: 
          variables.append(var.group(0)\
             .replace("VARIABLES = ","")\
             .replace("\"","")\
             .replace(" ","")\
             .split(','))
       except: return 0
   return array(variables).ravel()

def grid_data(df,resolution=0.1,variables=[]):
   from numpy import nan,linspace,meshgrid,arange,zeros
   from scipy.interpolate import griddata
   from pandas import DataFrame
   import os
   df_structured = DataFrame()

   if not len(variables):
      variables = df.columns

   xmin = df.x.min()
   xmax = df.x.max()
   ymin = df.y.min()
   ymax = df.y.max()

   # Create the structured grid that fits in this geometry
   xn = arange(xmin,xmax,resolution)
   yn = arange(ymin,ymax,resolution)
   #print "   Forming grid of {0} by {1}".format(len(xn),len(yn))
   X,Y = meshgrid(xn,yn)

   df_structured['x'] = X.ravel()
   df_structured['y'] = Y.ravel()

   # Fill this structured grid by interpolating values over the new nodes
   # for all the variables found in the above rotated data frame
   for v in variables:
      if v != 'x' and v != 'y':
         print "      Interpolating variable {0}".format(v)
         points = zeros((len(df['x'].values),2))
         values = zeros(len(df['x'].values))

         points[:,0] = df['x'].values
         points[:,1] = df['y'].values

         for i in range(len(df['x'].values)):
            values[i] = df.ix[
                  (df['x'] == points[i,0]) & 
                  (df['y'] == points[i,1]) ,
                  v]#.replace(nan,0)


         df_structured[v] = griddata(
               points,
               values,
               (df_structured['x'].values,df_structured['y'].values),
               method='linear'
               )
   return df_structured

def write_tecplot(case_name,df,outfile,interpolate=False):
   from numpy import nan,linspace,meshgrid,arange,zeros
   from scipy.interpolate import griddata
   from pandas import DataFrame
   import os

   variables = df.columns

   device,phi,alpha,U,loc = get_case_details(case_name)

    # Create the structured mesh for the new stitched mesh
   if interpolate:
      df = grid_data(df)
   # Construct the TECPLOT header for the combined data file
   header1 = "TITLE = \"{0}_phi{1}_alpha{2}_U{3}_loc{4}\"".format(device,phi,alpha,U,loc)
   header2 = "VARIABLES = "
   for var in df.columns:
           header2 = header2+"\"{0}\"".format(var)
           if var != df.columns[-1]:
                   header2 = header2+", "
   header3 = "ZONE T=\"Frame 1\" I={0}, J={1}, F=POINT\n".format(len(df['x'].unique()),len(df['y'].unique()))
   df = df.replace(nan,"NaN")
   df = df.sort(['y','x'])
   if os.path.isfile(outfile):
           os.remove(outfile)
   f = open(outfile,'w')
   f.write(header1+"\n")
   f.write(header2+"\n")
   f.write(header3+"\n")
   df.to_csv(f,sep=" ",header=False,index=False)
   f.close()


def concat_tecplot(case,out=False):
   """ Reads a series of TECPLOT files that come from a time averaged procedure from DaVis and 
   saves a pickle with a data frame in it for later processing

   Input: array of tecplot files
   """
   import os
   from numpy import array,unique,nan
   import pandas as pd
   import pickle


   tecplot_files = [f for f in os.listdir(case) if os.path.splitext(f)[1]=='.dat']

   # Get the uniquely different variables found in the tecplot files

   variables = []
   for tec in tecplot_files:
      var = get_variables_in_tecplot(os.path.join(case,tec))
      variables.append(var)

   variables = unique(array(variables).ravel())

   # Form the data frame

   for tec in tecplot_files:
      # Get the coordinate points
      var = get_variables_in_tecplot(os.path.join(case,tec))
      df_tmp = pd.read_table(
            os.path.join(case,tec),
            skiprows=3,
            sep=' ',
            names=var
            )
      if tec == tecplot_files[0]:
         df = pd.DataFrame()
         df['x'] = df_tmp['x']
         df['y'] = df_tmp['y']

      # Get the unique variable in this tecplot file
      u_var = [f for f in var if not f=='x' and not f=='y'][0]
      df[u_var] = df_tmp[u_var]

      # Get the streamwise vorticity info
   #df['Vorticity_xy'] = get_vorticity(df)

   try:
       device,phi,alpha,U,loc,reprocessed = get_case_details_from_filename(os.path.split(case)[-1])
   except:
       device,phi,alpha,U,loc,reprocessed = get_case_details_from_filename(case.split("/")[-2])

   if reprocessed:
      rep = "_Reprocessed"
   else:
      rep = ''
   if out:
       out_filename = out
   else:
       out_filename = "{0}_phi{1}_alpha{2}_U{3}_loc{4}{5}".format(device,phi,alpha,U,loc,rep)

   # Construct the TECPLOT header for the combined data file
   header1 = "TITLE = \"{0}_phi{1}_alpha{2}_U{3}_loc{4}\"".format(device,phi,alpha,U,loc)
   header2 = "VARIABLES = "
   for var in df.columns:
           header2 = header2+"\"{0}\"".format(var)
           if var != df.columns[-1]:
                   header2 = header2+", "
   fp = open(os.path.join(case,tecplot_files[0]))
   for i, line in enumerate(fp):
       if i == 2:
           header3 = line
   df = df.replace(nan,"NaN")
   if os.path.isfile(os.path.join('ConvertedData',out_filename+'.dat')):
           os.remove(os.path.join('ConvertedData',out_filename+'.dat'))
   if not os.path.isdir('ConvertedData'):
           os.makedirs('ConvertedData')
   f = open(os.path.join('ConvertedData',out_filename+'.dat'),'a')
   f.write(header1+"\n")
   f.write(header2+"\n")
   f.write(header3+"\n")
   df.to_csv(f,sep=" ",header=False,index=False)
   f.close()

   try:
      pickle.dump(
            df,
            open(out,'wb')
            )
      return 1
   except:
      return 0


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

def get_case_details(case_name):
   from re import findall
   device  = case_name.split("_")[0]
   phi = findall("phi-?[0-9]",case_name)
   if len(phi):
      phi = phi[0].replace("phi","")
   else:
      phi = "NA"
   alpha = findall("alpha-?[0-9][0-9]?",case_name)
   if len(alpha):
      alpha = alpha[0].replace("alpha","")
   else:
      alpha = "NA"
   U = findall("U-?[0-9]+",case_name)
   if len(U):
      U = U[0].replace("U","")
   else:
      U = "NA"
   loc = findall("loc[0-9][0-9]",case_name)
   if len(loc):
      loc = loc[0].replace("loc","")
   else:
      loc = "NA"
   return device,phi,alpha,U,loc

def get_vorticity(df):
    from numpy import meshgrid,shape,zeros
    from sys import exit

    if "Vorticity_xy" in df.columns:
        # Do nothing and return the same DF
        return df

    # Get shape of 2D plane
    nx = len(df['x'].unique())
    ny = len(df['y'].unique())
    Ux = df['Avg_Vx'].values.reshape((ny,nx))
    Uy = df['Avg_Vy'].values.reshape((ny,nx))
    ax = df['x'].values.reshape((ny,nx))/1000. # [mm] -> [m]
    ay = df['y'].values.reshape((ny,nx))/1000. # [mm] -> [m]

    i,j = shape(Ux)

    # Sanity check:
    if i != shape(Uy)[0] or i != shape(ax)[0] or i != shape(ay)[0]:
        exit("   The shape of the matrices while getting the vorticity is not the same!")

    duy_dax = zeros((i,j))
    dux_day = zeros((i,j))
    for ii in range(1,i-1):
        for jj in range(1,j-1):
            duy_dax[ii,jj] = (Uy[ii,jj+1]-Uy[ii,jj-1])/(ax[ii,jj+1]-ax[ii,jj-1])
    for ii in range(1,i-1):
        for jj in range(1,j-1):
            dux_day[ii,jj] = (Ux[ii+1,jj]-Ux[ii-1,jj])/(ay[ii+1,jj]-ay[ii-1,jj])

    vorticity = duy_dax - dux_day

    df['Vorticity_xy'] = vorticity.ravel()

    return df

def plot_airfoil(case,ax,factor=1):
   from numpy import array

   mask = rotate_mask(case)

   x = []; y = []
   if "STE" in case:
      x = [mask[0][0],mask[1][0],mask[-2][0],mask[-1][0]]
      y = [mask[0][1],mask[1][1],mask[-2][1],mask[-1][1]]
   else:
      for m in mask:
         x.append(m[0])
         y.append(m[1])
   ax.plot(array(x)*factor,array(y)*factor,lw=7,c='k')
   return ax

def get_mask_angle(case,to_freestream=True):
   import Masks as masks
   from math import atan,cos,sin,tan
   from numpy import deg2rad,zeros
   import Functions as piv
   from copy import copy

   device,phi,alpha,U,loc = piv.get_case_details(case)

   if to_freestream:
      incidence_correction = 0
   else:
      incidence_correction = abs(float(alpha))+abs(float(phi))

   angle = -atan( (masks.Masks[case][2][0] - masks.Masks[case][1][0]) \
         / (masks.Masks[case][2][1] - masks.Masks[case][1][1]) )\
         + deg2rad(90+incidence_correction)

   return angle

def rotate_mask(case):
   import Masks as masks
   from math import atan,cos,sin,tan
   from numpy import deg2rad,zeros
   import Functions as piv
   from copy import copy

   angle = get_mask_angle(case)
   device,phi,alpha,U,loc = piv.get_case_details(case)
   
   rotated_mask = zeros((6,2))
   for m,i in zip(masks.Masks[case],range(5)):
      cmx = m[0] - masks.Masks[case][1][0]
      cmy = m[1] - masks.Masks[case][1][1]
      x =  cmx*cos(angle) + cmy*sin(angle)
      y = -cmx*sin(angle) + cmy*cos(angle) 
      rotated_mask[i] = [x,y]
   rotated_mask[2] = [
         40.*sin(deg2rad(90-float(alpha)-float(phi))),
         -40.*cos(deg2rad(90-float(alpha)-float(phi)))
         ]
   rotated_mask[3] = [
         rotated_mask[2][0]-cos(deg2rad(90-float(alpha)-float(phi))),
         rotated_mask[2][1]-sin(deg2rad(90-float(alpha)-float(phi))),
         ]
   rotated_mask[4] = [
         rotated_mask[1][0]-cos(deg2rad(90-float(alpha)-float(phi))),
         rotated_mask[1][1]-sin(deg2rad(90-float(alpha)-float(phi))),
         ]
   rotated_mask[5][0] = rotated_mask[0][0]
   rotated_mask[5][1] = rotated_mask[4][1] - rotated_mask[5][0]*tan(deg2rad(float(alpha)-7)) 
   return rotated_mask


def make_surface_plot(cases=False,variable=False,output='test.png',
                      shape=False,zlabel='',case_names=False,
                      U=False,levels=False,subtract_mean = False,
                      ticks=[],streamlines=False,stream_factor=1,
                      root='/',figsize=False,fontsize=30):
    from matplotlib import pyplot as plt
    import matplotlib as mpl
    from numpy import meshgrid
    import re
    import os
    from matplotlib import rc

    rc('text',usetex=True)
    rc('font',family='serif')

    base_cmap = plt.get_cmap("RdYlBu_r")

    if not cases or not variable or not shape:
        return 0

    figsize_factor = 0.5
    if not figsize:
       if shape[1] == 2:
          figsize = (figsize_factor*12,figsize_factor*10)
       else:
          figsize = (figsize_factor*20,figsize_factor*10)
    fig,axes = plt.subplots(ncols=shape[1],nrows=shape[0], sharex=True, sharey=True,figsize=figsize)

    for case,ax,i in zip(cases,axes.flat,range(len(cases))):
        df = read_tecplot_file(os.path.join(root,case))
        case_name = re.findall(
            "[A-Za-z0-9]+_phi-?[0-9]_alpha-?[0-9]",
            os.path.splitext(case)[0])[0] 
        if subtract_mean: mean = df[variable].mean()
        else: mean = 0
        if variable == "Vorticity_xy":
           df = get_vorticity(df)
        if U:
            df["Avg_Vx"] = (df["Avg_Vx"]-mean)/float(U)
            df["Avg_Vy"] = (df["Avg_Vy"]-mean)/float(U)
            df["Avg_Vz"] = (df["Avg_Vz"]-mean)/float(U)
            df["RMS_Vx"] = (df["RMS_Vx"]-mean)/float(U)
            df["RMS_Vy"] = (df["RMS_Vy"]-mean)/float(U)
            df["RMS_Vz"] = (df["RMS_Vz"]-mean)/float(U)
        X,Y  = meshgrid(df.x.unique()/40., df.y.unique()/40.)
        try:
            Z = df[variable].reshape(X.shape)
        except KeyError:
           try:
              Z = df[DaVis_naming_convention_dictionary[variable]].reshape(X.shape)
           except KeyError:
              print "   ERROR: Variable {0} not available amongst\n\t{1}".format(variable,df.columns)
              raise
        ax.set_ylim(-.5,.5)
        ax.set_xlim(-10/40.,1.)
        ax.set(aspect=1, adjustable='box-forced')
        plot_airfoil(case,ax,factor=1/40.)
        contf = ax.contourf(X,Y,Z,cmap = base_cmap,levels=levels)
        if streamlines:
           col = 'w'
        else:
           col = 'k'
        cont  = ax.contour(X,Y,Z,colors=col,levels=levels)
        for col in cont.collections:
            col.set_linestyle('solid')
        if streamlines:
            ax.streamplot(
                  X, 
                  Y, 
                  df['Avg_Vx'].reshape(X.shape)*stream_factor, 
                  df['Avg_Vy'].reshape(X.shape),
                  color='k',linewidth=df['Avg_Vx'].values.reshape(X.shape)*3.,
                  density=[0.75,0.75]
                  )
        if case_names:
            ax.annotate(case_names[i], xy=(0.5, 1.03), xycoords='axes fraction', fontsize=fontsize,
                horizontalalignment='center', verticalalignment='bottom')

        plt.setp(ax.get_xticklabels(), rotation=45,fontsize=fontsize)
        plt.setp(ax.get_yticklabels(), fontsize=fontsize)
    plt.subplots_adjust(hspace=0.0)
    if shape[1] == 3:
        y_loc = 0.00
        x_loc = 0.08
        #y_loc = -0.08
        #x_loc = 0.0
    elif shape[1] == 4:
        y_loc = 0.04
        x_loc = 0.08
    elif shape[1] == 2:
        y_loc = 0.04
        x_loc = 0.08
        if shape[0] == 2:
           y_loc = 0.05
           x_loc = 0.08
    else:
        y_loc = 0.0
        x_loc = 0.0
    fig.text((1./shape[1])*(shape[1]/2.-0.05),y_loc,"$x/2h$",ha='center',fontsize=fontsize)
    fig.text(x_loc,0.5,"$y/2h$",ha='center',va='center',fontsize=fontsize,rotation=90)
    clb = plt.colorbar(contf, ax=axes.ravel().tolist(),fraction=5.*0.034/shape[1], pad=0.04)
    clb.set_label(zlabel,fontsize=fontsize,labelpad=30)
    cl = plt.getp(clb.ax, 'ymajorticklabels') 
    plt.setp(cl, fontsize=fontsize)
    if len(ticks):
       clb.set_ticks(ticks)
    fig.savefig(output, bbox_inches='tight')
    return 1


def rotate_tecplot_file(tecplot_file,angle,center):
    from copy import copy
    from math import sin,cos
    import pandas as pd
    from pandas import DataFrame
    from scipy.interpolate import interp2d,griddata
    from numpy import linspace,arange,meshgrid,zeros,rad2deg,nan,sqrt
    from os.path import split
    from sys import exit

    if type(tecplot_file) == pd.core.frame.DataFrame:
       df = tecplot_file
    else:
       df = read_tecplot_file(tecplot_file)

    # Apply mask
    #df = apply_mask(df,split(tecplot_file)[-1])

    x = df['x'] - [center[0]]*len(df['x'].values)
    y = df['y'] - [center[1]]*len(df['y'].values)
    df_rotated = DataFrame()
    print "   Rotating by {0}".format(rad2deg(angle))
    df_rotated['x'] =  x*cos(angle) + y*sin(angle) 
    df_rotated['y'] = -x*sin(angle) + y*cos(angle) 

    # Velocity components
    df_rotated['Avg_Vx'] =  df['Avg_Vx']*cos(angle) + df['Avg_Vy']*sin(angle)
    df_rotated['Avg_Vy'] = -df['Avg_Vx']*sin(angle) + df['Avg_Vy']*cos(angle)
    df_rotated['Avg_Vz'] = df['Avg_Vz']

    # RMS components
    try:
       df_rotated['Reynold_stress_XX'] = df['Reynold_stress_XX']*cos(angle)**2\
               + df['Reynold_stress_YY']*sin(angle)**2\
               + 2.*cos(angle)*sin(angle)*df['Reynold_stress_XY']

       df_rotated['Reynold_stress_YY'] = df['Reynold_stress_XX']*sin(angle)**2\
               + df['Reynold_stress_YY']*cos(angle)**2\
               - 2.*cos(angle)*sin(angle)*df['Reynold_stress_XY']
       df_rotated['Reynold_stress_XY'] = df['Reynold_stress_XY']*cos(2*angle) \
             + (cos(angle)*sin(angle)) * ( df['Reynold_stress_XX']**2-df['Reynold_stress_YY']**2)

       df_rotated['RMS_Vx'] = sqrt(df_rotated['Reynold_stress_XX'].replace(nan,0).values)
       df_rotated['RMS_Vy'] = sqrt(df_rotated['Reynold_stress_YY'].replace(nan,0).values)
       df_rotated['RMS_Vx'] = df_rotated['RMS_Vx'].replace(0,nan)
       df_rotated['RMS_Vy'] = df_rotated['RMS_Vy'].replace(0,nan)
    except KeyError:
       print df.columns; raise

    try:
        df_rotated['RMS_Vz'] = df['RMSVz']
    except KeyError:
        df_rotated['RMS_Vz'] = df['RMS_Vz']

    return df_rotated

def stitch_cases(caseLS,caseRS):

    from pandas import DataFrame
    from numpy import linspace,meshgrid,zeros,nan,rad2deg,deg2rad
    import Masks as masks
    from math import atan
    import os


    case_name_LS = os.path.split(caseLS)[-1]
    case_name_RS = os.path.split(caseRS)[-1]
    device,phi,alpha,U,loc = get_case_details(case_name_RS)
    # Get the angles so that the frames are rotated in respect to the serration
    # surface to make it vertical
    angleLS = -atan( (masks.Masks[case_name_LS][2][0] - masks.Masks[case_name_LS][1][0]) \
          / (masks.Masks[case_name_LS][2][1] - masks.Masks[case_name_LS][1][1]) )
    #print "   For points ({0:.2f},{1:.2f}) and ({2:.2f},{3:.2f}), rotating by {4:.2f}".format(
    #      masks.Masks[case_name_LS][2][0] , masks.Masks[case_name_LS][2][1] ,
    #      masks.Masks[case_name_LS][1][0] , masks.Masks[case_name_LS][1][1] ,
    #      rad2deg(angleLS)
    #      )


    angleRS = -atan( (masks.Masks[case_name_RS][2][0] - masks.Masks[case_name_RS][1][0]) \
          / (masks.Masks[case_name_RS][2][1] - masks.Masks[case_name_RS][1][1]) )
    #print "   For points ({0:.2f},{1:.2f}) and ({2:.2f},{3:.2f}), rotating by {4:.2f}".format(
    #      masks.Masks[case_name_RS][2][0] , masks.Masks[case_name_RS][2][1] ,
    #      masks.Masks[case_name_RS][1][0] , masks.Masks[case_name_RS][1][1] ,
    #      rad2deg(angleRS)
    #      )

    # Apply the masks to the data frames
    caseLS_df = rotate_tecplot_file(
          apply_mask(
             read_tecplot_file(caseLS),
             case_name_LS
             ),
          angleLS,
          masks.Masks[case_name_LS][1]
          )
    caseRS_df = rotate_tecplot_file(
          apply_mask(
             read_tecplot_file(caseRS),
             case_name_RS
             ),
          angleRS,
          masks.Masks[case_name_RS][1]
          )

    # Center the cases in the [1] geometry position
    #mask = masks.Masks[os.path.split(caseLS)[-1]]
    #caseLS_df['x'] = caseLS_df['x'] - [mask[1][0]]*len(caseLS_df['x'].values)
    #caseLS_df['y'] = caseLS_df['y'] - [mask[1][1]]*len(caseLS_df['y'].values)
    #mask = masks.Masks[os.path.split(caseRS)[-1]]
    #caseRS_df['x'] = caseRS_df['x'] - [mask[1][0]]*len(caseRS_df['x'].values)
    #caseRS_df['y'] = caseRS_df['y'] - [mask[1][1]]*len(caseRS_df['y'].values)

    # Only grab the parts "above" the serration line
    df_ls = caseLS_df[caseLS_df['x']<0]
    df_rs = caseRS_df[caseRS_df['x']<0]

    mirror_variables = ["x","Avg_Vx"]
    for c in caseRS_df.columns:
        if c in mirror_variables:
           try:
               val = df_ls[c].values*([-1]*len(df_ls[c].values))
               df_ls.loc[:,c] = val
           except KeyError:
               print "   Error: looking for variable {0} but couldn't find it in {1}".format(mv,df_ls.columns)

    #return df_ls.append(df_rs[df_rs['y']<0].sort(['y','x'],ascending=False),ignore_index=True)
    stitched_df = df_rs.append(df_ls.sort(['y','x'],ascending=False),ignore_index=True)

    stitched_df = rotate_tecplot_file(
          stitched_df,
          deg2rad(float(alpha)+float(phi)+90),
          [0,0]
          )

    return stitched_df.replace(0,nan),df_ls.replace(0,nan),df_rs.replace(0,nan)

def read_tecplot_file(tecplot_file):
    """Reads in a tecplot file, given, and returns a pandas data frame

    """
    import pandas as pd 
    import re

    # Get available variables
    f = open(tecplot_file,'ro')
    var_flag = False
    dev_flag = False
    for line in f:
        string = re.findall("^VARIABLES[ _A-Za-z0-9,\"=]+",line)
        if string:
            variables = [v.replace(' ','_').replace("\"","") \
                         for v in string[0].replace("VARIABLES = ",'')\
                         .split(", ")]
            variables = [v for v in variables if len(v)]
            var_flag = True
        string = re.findall("^TITLE = [ -_A-Za-z0-9,\"=]+",line)
        if string:
            dev_flag = True
        if var_flag and dev_flag:
            break
    f.close()

    for v in range(len(variables)):
        if variables[v] in DaVis_naming_convention_dictionary.keys():
            variables[v] = DaVis_naming_convention_dictionary[variables[v]]\
                    .replace(" ","_")

    # Correct for error in the "Flow angle" variable name in the tecplot files
    variables = [v.replace("Flow_angle","") for v in variables]
    variables.append("Flow_angle")

    # Put the data into a data frame
    data = pd.read_table(
            tecplot_file,
            skiprows  = 4,
            names     = variables,
            sep       = '[ \t]+',
            index_col = False
            )

    return data
