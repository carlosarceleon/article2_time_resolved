from pandas import read_pickle
from article2_data_analysis_routines import return_bl_parameters

def quad_func( x, a , b, c):
    return a*x**2 + b*x + c

def plot_all_bls( bl_pickle_file = 0 ):
    import matplotlib.pyplot as plt
    from matplotlib import rc
    import seaborn as sns
    from pandas import read_pickle
    from scipy.optimize import curve_fit

    rc('text',usetex=True)
    rc('font',weight='normal')

    sns.set_context('paper')
    sns.set(font='serif',font_scale=1.6,style='whitegrid',
            rc={"axes.axisbelow": False})
    rc('font',family='serif', serif='cm10')

    if not bl_pickle_file:
        bl_pickle_file = 'BLData_staged.p'

    bl_df = read_pickle( bl_pickle_file )

    alphas = bl_df.alpha.unique()

    for a in alphas:
        fig,axes = plt.subplots( 4, 1, sharex = True , figsize = (15,12))
        cases = bl_df[ bl_df.alpha == a ].case.unique()

        cnt = 0
        variables = [ 'Ue', 'delta_99', 'delta_displacement', 
                    'delta_momentum' ]
        for var in variables:
            for c in cases:

                if 'z00' in c and not "STE" in c:
                    limits = (0, 40)
                elif 'z05' in c:
                    limits = (0, 20)
                else:
                    limits = (-8, 10)

                case_df = bl_df[ 
                    ( bl_df.case == c ) & \
                    ( bl_df.alpha == a )
                ][ ['x',var] ].dropna()

                skip = False
                try:
                    popt, pvar = curve_fit( 
                        quad_func, 
                        case_df[ (case_df.x < limits[1] ) & \
                                ( case_df.x > limits[0] )].x, 
                        case_df[ (case_df.x < limits[1] ) & \
                                ( case_df.x > limits[0] )][var])
                except TypeError:
                    print c
                    skip = True

                axes[cnt].plot( 
                    case_df[ case_df.x < limits[1] ].x,
                    case_df[ case_df.x < limits[1] ][var],
                    'o',
                    label = c.replace("_"," ")
                )

                if not skip:
                    axes[cnt].plot( 
                        case_df[ case_df.x < limits[1] ].x,
                        quad_func( case_df[ case_df.x < limits[1] ].x, 
                                  popt[0], popt[1], popt[2] ),
                        '-',
                        c = 'k',
                        lw = 2
                    )

            cnt += 1

        for ax,var in zip(axes,variables):
            ax.set_xlim( -5, 40 )
            ax.set_xlabel('$x$ [mm]')
            ax.set_ylabel(var.replace('_',' ')+ ' [mm]')

        axes[0].legend( bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
           ncol=2, mode="expand", borderaxespad=0.)

        plt.savefig("BL_a{0}".format(a), bbox_inches = 'tight')

def decript_case_name(case_name):
    from re import findall

    try:
        trailing_edge = findall('[litSrTE201R]+', case_name)[0]
    except IndexError:
        print case_name

    phi = 0 ; alpha = 0 ; U = 0 ; z = 0

    try:
        phi   = findall('phi[0-9]',case_name)[0].replace('phi','')
        alpha = findall('alpha-?[0-9][0-9]?',case_name)[0].replace('alpha','')
        z     = findall('loc[0-9][0-9]',case_name)[0].replace('loc','')

        if   z == '00' : z = '0'
        elif z == '05' : z = '0.25'
        elif z == '10' : z = '0.5'

    except:
        # In case it's the DaVis naming convention ############################
        if not phi:
            phi   = findall('p-?[0-9]',case_name)[0].replace('p','')
        if not alpha:
            alpha = findall('a-?[0-9][0-9]?',case_name)[0].replace('a','')
        if not z:
            z     = findall('z[0-9][0-9]',case_name)[0].replace('z','')
        # #####################################################################

    U = findall('U[0-9][0-9]',case_name)[0].replace('U','')

    return trailing_edge,phi,alpha,U,z

def clean_data( df , var , window = 3, threshold = 0.5):
    from pandas import rolling_median
    from numpy import abs

    original_columns = df.columns

    rolling_median = rolling_median(
        df[var],
        window = window,
        center = True,
    )

    df['this_mean'] = ( df[var].mean() + rolling_median ) / 2.

    df = df[ abs( df[var] - df.this_mean ) <= threshold * df[var].std() ]

    return df[ original_columns ]



def get_edge_velocity(df, 
                      condition = 'vorticity_integration_rate_of_change', 
                      threshold = 0.01,
                      Ue_file = 'Ue_data.csv',
                     ):
    from scipy.integrate import simps
    from numpy import diff,nan
    
    U_edge = 0 

    df = df.sort_values( by = 'y', ascending = True )
    df = df[df.y>5]
    df = df.reset_index(drop=True)

    # Get the integrated vorticity #############################################
    integration_result = []
    for i in range(len(df.y)):
        integration_result.append(
            - simps(df.vorticity_xy.ix[:i], 
                  x=df.y.ix[:i]/1000., 
                  even='avg'
                 )
        )
    if not len( integration_result ):
        return nan

    df['vorticity_integration'] = integration_result
    rate_of_change = diff(df.vorticity_integration)/diff(df.y/1000.)

    if not len( rate_of_change ):
        return nan

    rate_of_change = list(rate_of_change) + [rate_of_change[-1]]
    df['vorticity_integration_rate_of_change'] = rate_of_change
            
    ############################################################################

    # Get the rms rates of vertical change #####################################
    df['u_rms_rate_of_change'] = list(diff(df.u_rms)/diff(df.y/1000.))+[0]
    df['v_rms_rate_of_change'] = list(diff(df.v_rms)/diff(df.y/1000.))+[0]
    ############################################################################

    #df.plot( condition, 'y' )
    #plt.show()

    while not U_edge:
        for y_loc in df.iterrows():
            if y_loc[1][condition] < threshold:
                U_edge = y_loc[1]['u']
                break

        if not U_edge:
            threshold += 0.01
            if threshold > 2:
                print "      Could not find a U_e for x = {0:.2f}, {1}".format(
                    df.x.unique()[0], df.case_name.unique()[0]
                )
                return nan

    return U_edge


def get_boundary_layer_values( df , U_e = 0):
    import pandas as pd

    def get_delta_99(df,U_e):
        from numpy import isnan,nan

        if isnan(U_e):
            return nan

        df = df.append( {'u' : 0.99*U_e} , ignore_index = True)
        df = df.sort_values( by = 'u' )
        df = df.apply(pd.Series.interpolate)

        for y in sorted(df.y.unique()):
            if df[ df.y == y ].u.values[0] > 0.99*U_e:
                break

        #return df[ df.u == 0.99*U_e ].y.values[0]
        return y

    def get_delta_momentum(df,U_e, delta_99):
        from scipy.integrate import simps
        from numpy import isnan,nan

        df = df[ ( df.y > 1.0 ) & ( df.y < delta_99 ) ]

        if isnan(U_e):
            return nan
        return simps(
            ( df.u / U_e ) * ( 1 - df.u / U_e ), 
            x = df.y , 
            even='avg'
        )

    def get_delta_displacement(df,U_e, delta_99):
        from scipy.integrate import simps
        from numpy import isnan,nan

        df = df[ ( df.y > 1.0 ) & ( df.y < delta_99 ) ]
        if df.empty:
            return 0

        if isnan(U_e):
            return nan
        return simps(
            1 - df.u / U_e , 
            x = df.y , 
            even = 'avg'
        )

    if not U_e:
        U_e = get_edge_velocity(df)

    delta_99           = get_delta_99(df,U_e)
    delta_displacement = get_delta_displacement(df,U_e, delta_99)
    delta_momentum     = get_delta_momentum(df,U_e, delta_99)

    return U_e, delta_99, delta_displacement, delta_momentum


def run_bl_analysis( pickles_folder = 0 ):
    import matplotlib.pyplot as plt
    from os                           import listdir
    from os.path                      import join
    from pandas                       import read_pickle, DataFrame
    from article2_time_resolved_routines import find_nearest

    if not pickles_folder:
        pickles_folder = '/home/carlos/Documents/PhD/Articles/'+\
                'Article_3/Scripts/time_resolved/averaged_data'

    case_pickles = [
        f for f in listdir( pickles_folder ) if f.endswith(".p") \
        if not 'Slit' in f and 'alpha0' in f and 'phi0' in f\
        and not "mean_flow_rotated" in f
    ]

    bl_df = DataFrame()

    fig, ax = plt.subplots( 1, 1 )

    for cp in case_pickles:
        case_bl_df = DataFrame()

        df = read_pickle( join( pickles_folder, cp ) )
        df = df.sort_values( by = [ 'x', 'y' ] )

        if 'loc00' in cp and not 'STE' in cp:
            x_bl_loc = 38
        elif 'loc05' in cp:
            x_bl_loc = 18
        elif 'loc10' in cp or 'STE' in cp:
            x_bl_loc = -2

        available_x_loc = find_nearest( x_bl_loc, df.x.values )

        trailing_edge,phi,alpha,U,z = \
                decript_case_name(cp)

        case_name = "{0}_a{1}_p{2}_U20_z{3:02.0f}_tr".\
                format( trailing_edge, alpha, phi, float(z)*20 )

        print "   Running {0}".format(case_name)

        # First get the edge velocity, because it needs to be cleaned up a bit #
        ue_df = DataFrame()
        for x in df.x.unique():
            local_x_df = df[ ( df.x == x ) & ( df.y >= 0 ) ]

            ue_df = ue_df.append(
                { 'U_e' : get_edge_velocity( local_x_df ),
                 'x' : x}, ignore_index = True
            )
        # ######################################################################

        ue_df = clean_data( ue_df, 'U_e' , window = 10, threshold = 1.0 )

        for x , U_e_x in zip( ue_df.x.values, ue_df.U_e.values ):
            local_x_df = df[ ( df.x == x ) & ( df.y >= 0 ) & ( df.y < 20 ) ]

            if x == available_x_loc:
                ax.plot( 
                    local_x_df.u / U_e_x,
                    local_x_df.y, 
                    label = cp.replace("_"," ") 
                )

            U_e_loc, delta_99, delta_displacement, delta_momentum = \
                    get_boundary_layer_values( local_x_df, U_e_x )

            data = {
                'case':               case_name,
                'Ue':                 U_e_x,
                'delta_99':           delta_99,
                'delta_displacement': delta_displacement,
                'delta_momentum':     delta_momentum,
                'x':                  x,
                'trailing_edge':      trailing_edge,
                'phi':                phi,
                'alpha':              alpha,
                'z':                  z
            }

            case_bl_df = case_bl_df.append( 
                DataFrame( data, index = [0] ),
                ignore_index = True
            )

        if 'delta_99' in case_bl_df.columns:
            case_bl_df = clean_data( case_bl_df, 'delta_99', window = 10 , 
                                    threshold = 1.0 )

        bl_df = bl_df.append( case_bl_df, ignore_index = True )

    bl_df.to_pickle("BLData_staged.p")

    plt.legend( loc = 'best' )
    plt.xlim( 0, 1 )
    plt.savefig( "InterestingBLs.png" )


bldf = read_pickle( 'BLData_staged.p' )


cases = [
    'STE_a0_p0_U20_z00_tr', 
    'Sr20R21_a0_p0_U20_z10_tr', 
    'Sr20R21_a0_p0_U20_z05_tr', 
    'Sr20R21_a0_p0_U20_z00_tr'
]

x_locs = [ 0, 0, 20, 40 ]

for c,x in zip( cases, x_locs ):
    print c
    print return_bl_parameters( c, [x] , bl_pickle_file = 'BLData_staged.p' )

#plot_all_bls( 'BLData.p' )
