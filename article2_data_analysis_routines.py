def get_bl_parameters(
    df, 
    bl_file = 'Boundary_layer_information.csv'
):
    import pandas as pd
    bl_data = pd.read_csv( bl_file )

    if 'file' in df.columns:
        case_col_name = 'file'
    elif 'case' in df.columns:
        case_col_name = 'case'
    else:
        print df
    # Device ###################################################################
    if 'Sr20R21' in df[case_col_name].unique()[0]:
        device = 'serrated'
    elif "STE" in df[case_col_name].unique()[0]:
        device = 'straight'
    # Z location ###############################################################
    if 'z00' in df[case_col_name].unique()[0]:
        z_loc = 0
    elif 'z05' in df[case_col_name].unique()[0]:
        z_loc = 0.25
    elif 'z10' in df[case_col_name].unique()[0]:
        z_loc = 0.5
    # X location ###############################################################
    if 'near_x' in df.columns:
        x_loc = df.near_x.unique()[0]
    elif 'near_x_upwind' in df.columns:
        x_loc = df.near_x_upwind.unique()[0]

    bl_case_data = bl_data[
        ( bl_data.z_loc         == z_loc   ) &\
        ( bl_data.x_loc         >= x_loc-2 ) &\
        ( bl_data.x_loc         <= x_loc+2 ) &\
        ( bl_data.Trailing_edge == device  )
    ]

    if bl_case_data.empty:
        print df

    return bl_case_data



def remove_angle_jumps(df):
    from numpy import sign
    from math import pi

    for var in ['u','v']:
        df['phi_'+var].loc[df['phi_'+var]<0] = \
                df['phi_'+var].loc[df['phi_'+var]<0] + pi

        for ix in range(len(df))[:-2]:
            dif = df['phi_'+var].ix[ix+1] - df['phi_'+var].ix[ix]
            if abs(dif) > pi*0.4:
                df['phi_'+var].ix[ix+1] = df['phi_'+var].ix[ix+1] \
                        - sign(dif) * pi

        df['phi_'+var].loc[df['phi_'+var]<0] = \
                df['phi_'+var].loc[df['phi_'+var]<0] + pi

        #df['phi_'+var].loc[df['f_'+var] == df['f_'+var].max()] = \
        #        df['phi_'+var].loc[df['f_'+var] == df['f_'+var].max()] + pi

    return df

def calculate_Uc(df,delta_x):
    from math import pi
    from scipy.stats import linregress

    for var in ['u','v']:
        df = df.sort_values(by = [ 'f_' + var ], ascending=True )\
                .reset_index( drop = True )

        r_value  = 0
        consider = len(df)
        while r_value**2<0.97:
            df = df.sort_values(by = ['f_'+var]).\
                    ix[:consider].\
                    reset_index(drop=True)

            slope, intercept, r_value, p_value, std_err = linregress(
                df['f_'+var],
                df['phi_'+var],
            )
            consider -= 1
        
        if var == 'u':
            Uc_u = 2*pi*slope**-1*delta_x/1000.
            slope_u = slope
        if var == 'v':
            Uc_v = 2*pi*slope**-1*delta_x/1000.
            slope_v = slope

    return Uc_u, Uc_v, df, intercept, slope_u, slope_v

def do_the_coherence_analysis(df, schematic = ''):
    import pandas as pd

    coherence_df = pd.DataFrame()
    # Check how many cases there are ###########################################
    cases = df.case_name.unique()
    if not len(cases) > 1: 
        print " Found less than one case submitted"
        print "   {0}".cases
        return pd.DataFrame()

    for c in cases:

        # Split the upwind and downwind signals ################################
        x_locs = df[ df.case_name == c ].near_x.unique()

        if not len(x_locs) == 2: 
            print "   Didn't find more than two streamwise locations for"
            print "      {0}".format(c)
            print "   at y = {0}".format(df.near_y.unique())
            return pd.DataFrame()

        left_signal  = df[ (df.case_name == c) & (df.near_x == x_locs.min()) ]
        right_signal = df[ (df.case_name == c) & (df.near_x == x_locs.max()) ]
        # ######################################################################

        coherence_df = coherence_df.append(
            get_Uc_phi_and_coherence( left_signal, right_signal ),
            ignore_index = True
        )

    return coherence_df


def get_Uc_phi_and_coherence(signal1_df, signal2_df):
    import pandas as pd
    from scipy.signal import csd
    from numpy import abs,arctan,sqrt

    max_lag          = 10000
    freq_lower_limit = 300
    nperseg          = 2**6
    fs               = 10000
    x_1              = signal1_df.x.unique()[0]
    x_2              = signal2_df.x.unique()[0]
    y_1              = signal1_df.y.unique()[0]
    y_2              = signal2_df.y.unique()[0]

    delta_x          = sqrt(abs(x_1 - x_2)**2 + abs(y_1 - y_2)**2)

    df = pd.DataFrame()

    for var in ['u', 'v']:
        s1 = signal1_df[var].values[0:max_lag] \
                - signal1_df[var].values[0:max_lag].mean()
        s2 = signal2_df[var].values[0:max_lag] \
                - signal2_df[var].values[0:max_lag].mean()

        f,Pxy = csd(
            s2,s1,
            nperseg = nperseg,
            fs      = fs,
        )

        f,Pxx = csd(
            s1,s1,
            nperseg = nperseg,
            fs      = fs,
        )

        f,Pyy = csd(
            s2,s2,
            nperseg = nperseg,
            fs      = fs,
        )

        gamma_squared = abs( Pxy )**2 / ( Pxx * Pyy )

        gamma = sqrt(gamma_squared)

        Phi = arctan( Pxy.imag / Pxy.real )

        data = pd.DataFrame()
        data['f_'+var]     = f[:-1]
        data['phi_'+var]   = Phi[:-1]
        data['gamma_'+var] = gamma[:-1]
        data['mean_'+var]  = signal2_df[var].values[0:max_lag].mean()
        data['std_'+var]   = signal2_df[var].values[0:max_lag].std()


        data = data[
            data['f_'+var] >= freq_lower_limit
        ].reset_index( drop = True )

        df = pd.concat( [ df, data ], axis = 1 )

    df = remove_angle_jumps(df)

    Uc_u, Uc_v,data,intercept,slope_u,slope_v = calculate_Uc(
        df,
        delta_x = delta_x
    )

    df['Uc_u']            = Uc_u
    df['Uc_v']            = Uc_v
    df['trusted_f_v']     = data.f_v
    df['trusted_f_u']     = data.f_u
    df['slope_u']         = slope_u
    df['slope_v']         = slope_v
    df['x_upwind']        = signal1_df.x.unique()[0]
    df['near_x_upwind']   = signal1_df.near_x.unique()[0]
    df['y_upwind']        = signal1_df.y.unique()[0]
    df['near_y_upwind']   = signal1_df.near_y.unique()[0]
    df['x_downwind']      = signal2_df.x.unique()[0]
    df['near_x_downwind'] = signal2_df.near_x.unique()[0]
    df['y_downwind']      = signal2_df.y.unique()[0]
    df['near_y_downwind'] = signal2_df.near_y.unique()[0]
    df['case']            = signal1_df.case_name.unique()[0]

    return df

def get_color_and_marker(case_name):

    if "STE" in case_name:
        color  = (0.0, 0.4470588235294118, 0.6980392156862745)
        marker = 'x'
        cmap = 'Blues'
    elif 'z00' in case_name:
        color = (0.0, 0.6196078431372549, 0.45098039215686275)
        marker = '2'
        cmap = 'Greens'
    elif 'z05' in case_name:
        color = (0.8352941176470589, 0.3686274509803922, 0.0)
        marker = '+'
        cmap = 'Oranges'
    elif 'z10' in case_name:
        color = (0.8, 0.4745098039215686, 0.6549019607843137)
        marker = 'o'
        cmap = 'RdPu'

    else: print case_name; return 0,0
    return color,marker,cmap
                        
def do_the_reynolds_stress_quadrant_analysis(cases_df,y, plot_name = ''):

    from matplotlib import rc
    import seaborn as sns
    import matplotlib.pyplot as plt
    
    rc('text',usetex=True)
    rc('font',weight='normal')

    sns.set_context('paper')
    sns.set(font='serif',font_scale=3.0,style='whitegrid')
    rc('font',family='serif', serif='Linux Libertine')

    cases = sorted(cases_df.file.unique())
    fig,axes = plt.subplots(
        1,len(cases), 
        figsize = (figsize[0]*len(cases), figsize[1]), 
        sharex=True,
        sharey=True, 
    )
        
    for case_file, case_i, ax in zip(cases, range(len(cases)), axes):

        case = cases_df[cases_df.file == case_file]\
                .sort_values( by = ['time_step'] )

        bl_data = get_bl_parameters( case )

        color, marker, cmap = get_color_and_marker( 
            case_file
        )

        if not color and not marker:
            break

        case["uprime"] = ( case.u - case.u.mean() ) / bl_data.Ue.values[0]
        case["vprime"] = ( case.v - case.v.mean() ) / bl_data.Ue.values[0]
        
        ax.plot(
            case['uprime'].values[::50],
            case['vprime'].values[::50],
            ls              = '',
            marker          = marker,
            markeredgewidth = markeredgewidth,
            markerfacecolor = markerfacecolor,
            markeredgecolor = 'k',
            markersize      = 3*markersize/4.,
            mew             = mew,
            color           = 'k',
            alpha           = 0.3
        )

        kde = sns.kdeplot(case.uprime, case.vprime,
                    cmap         = cmap,
                    ax           = ax,
                    shade        = True,
                    shade_lowers = False,
                    gridsize     = 30
                   )
        
        ax.set_aspect('equal')
        kde.collections[0].set_alpha(0)

        ax.set_xlabel('')
        ax.set_ylabel('')

        ax.set_xlim( -0.3 , 0.3 )
        ax.set_ylim( -0.3 , 0.3 )

        ax.axhline( 0, ls = '--', lw=3 , c = 'k')
        ax.axvline( 0, ls = '--', lw=3 , c = 'k')

        ax.set_xlabel(r"$u'/u_e$")
        ax.grid(False)

    axes[0].set_ylabel(r"$v'/u_e$")

    axes[-1].text( 0.02, 0.2, r'$y/\delta_{{99}} = {0}$'.format(y/delta/1000.))


    plot_name = 'Results/ReynoldsQuadrant_{0}_ydelta{1:.2f}.png'\
        .format( plot_name,y/delta/1000. )

    fig.savefig( 
        plot_name.replace('.','_').replace('_png','.png'),
        bbox_inches = 'tight'
    )
    plt.cla()


def do_the_frequency_analysis(cases_df, y, plot_name = '', schematic = ''):
    
    from scipy.signal import welch
    from numpy import log10
    from matplotlib import rc
    import seaborn as sns
    import matplotlib.pyplot as plt
    import article2_time_resolved_routines as trr
    from matplotlib.cbook import get_sample_data
    
    rc('text',usetex=True)
    rc('font',weight='normal')

    sns.set_context('paper')
    sns.set(font='serif',font_scale=3.0,style='whitegrid')
    rc('font',family='serif', serif='Linux Libertine')

    if not len(cases_df):
        print "   No cases were passed to process"
        return 0

    for var in ['u','v']:

        figsize = (8,5)
        fig,ax = plt.subplots(1,1, figsize = figsize)
        
        kolmogorov_law_curve = [[0],[0]]

        if len(cases_df.file.unique()) <= 1:
            return 0

        for case_file in cases_df.file.unique():

            case = cases_df[cases_df.file == case_file]\
                    .sort_values( by = ['time_step'] )

            freq, Pxx = welch(
                x       = case[var].values,
                nperseg = nperseg,
                fs      = fs,
                scaling = 'spectrum',
            )

            case['case'] = case.file.unique()[0]
            bl_data = get_bl_parameters( case )

            color, marker, cmap = get_color_and_marker( 
                case_file
            )

            if not color and not marker:
                break

            res = dict( zip( 
                trr.get_Strouhal( freq, bl_data.delta_99.values[0], 
                                 bl_data.Ue.values[0] ),
                10 * log10( Pxx ),
            ))
            
            ax.plot(
                trr.get_Strouhal( freq, bl_data.delta_99.values[0], 
                                 bl_data.Ue.values[0] ),
                10 * log10( Pxx ),
                marker          = marker,
                markeredgewidth = markeredgewidth,
                markerfacecolor = markerfacecolor,
                markeredgecolor = color,
                markersize      = markersize,
                mew             = mew,
                color           = color
            )

            k_lims = ( sorted(res.keys())[8], sorted(res.keys())[20] )

            if not any(kolmogorov_law_curve[1]) or \
               res[k_lims[0]] > kolmogorov_law_curve[1][0]:

                kolmogorov_law_curve = get_kolmogorov_law_curve(x_lim = k_lims)
                kolmogorov_law_curve[1] = kolmogorov_law_curve[1] \
                        + res[k_lims[0]] + 3

        ax.plot( 
            kolmogorov_law_curve[0] , 
            kolmogorov_law_curve[1] , 
            '--',
            color = 'k' ,
            lw    = 3,
        )

        if not k_lims == (0,0):
            ax.text(
                k_lims[0]+(k_lims[1]-k_lims[0])/3.,
                kolmogorov_law_curve[1][0] - 3,
                "$\\textrm{St}^{-5/3}_\\delta$",
            )

        ax.set_xscale('log')
        

        ax.set_xlim( 0.09 , 2.2 )
        ax.set_ylim( -30 , 5 )
        if y//delta/1000. == 0.1:
            ax.set_xlabel(r"$\textrm{St_{99}} = f\delta_{99}/u_e$")
        ax.set_ylabel(
            r"$10\log_{10}\left(\Phi_{"+\
            str(var)+r"}\right)$ [dB]"
        )
        ax.text(
            1.8, -28,
            r'$y/\delta_{{99}} = {0}$'.format(y/delta/1000.),
            ha = 'right'
        )

        plt.grid(True, which='both')

        if schematic and y//delta/1000. == 0.1:
            im = plt.imread( get_sample_data( schematic  ) )
            newax = fig.add_axes([0.175, 0.175, 0.3, 0.3], anchor = 'SW', 
                                         zorder=100)
            newax.imshow(im)
            newax.axis('off')

        plot_composed_name = 'Results/FreqSpectra_{0}_ydelta{1:.2f}_{2}.png'\
            .format( plot_name, y/delta/1000., var )
        print " Going to save\n   {0}".format( plot_composed_name )

        fig.savefig( 
            plot_composed_name.replace('.','_').replace('_png','.png'),
            bbox_inches = 'tight'
        )
        plt.cla()

def get_kolmogorov_law_curve( x_lim = (0.5,1.5) ):
    from numpy import linspace, log10

    slope_x = linspace( x_lim[0], x_lim[1] ,100 )

    slope_y = 10*log10(slope_x**(-5/3.))

    return [slope_x,slope_y]

def plot_coherence_Uc_phi( coherence_df , plot_name = '', schematic = ''):
    import matplotlib.pyplot as plt
    from article2_time_resolved_routines import get_Strouhal
    from matplotlib import rc
    import seaborn as sns
    from numpy import array,linspace,exp
    from math import pi
    from matplotlib.cbook import get_sample_data


    rc('text',usetex=True)
    rc('font',weight='normal')

    sns.set_context('paper')
    sns.set(font='serif',font_scale=3.0,style='whitegrid')
    rc('font',family='serif', serif='Linux Libertine')

    for var in ['u','v']:

        fig_Uc,    ax_Uc    = plt.subplots(1, 1, figsize = figsize)
        fig_std,  ax_std  = plt.subplots(1, 1, figsize = figsize)

        for case in coherence_df.case.unique():

            case_df = coherence_df[coherence_df.case == case]
            
            bl_data = get_bl_parameters(case_df)

            Uc = []
            Um = []
            std = []
            y_locs = sorted(case_df\
                            .near_y_downwind.unique())

            real_y_locs = []
            for y_loc in y_locs:
                y_df = case_df[ 
                    (case_df.near_y_downwind == y_loc) &\
                    (case_df.case == case) 
                ]

                Uc.append(
                    y_df['Uc_'+var].unique()[0]
                )
                Um.append(
                    y_df['mean_u'].unique()[0]
                )
                std.append(
                    y_df['std_u'].unique()[0]
                )

                real_y_locs.append(
                    case_df[ case_df.near_y_downwind == y_loc ]\
                    .y_downwind.unique()[0]
                )

            real_y_locs = array(real_y_locs)

            color, marker, cmap = get_color_and_marker( case )

            plot_config = {
                'marker'          : marker,
                'markeredgewidth' : markeredgewidth,
                'markerfacecolor' : markerfacecolor,
                'markeredgecolor' : color,
                'markersize'      : markersize,
                'mew'             : mew,
                'color'           : color,
            }

            ax_std.plot(
                array(std) / bl_data.Ue.values[0], 
                real_y_locs / bl_data.delta_99.values[0],
                ls='',
                **plot_config
            )

            ax_Uc.plot(
                array(Uc) / bl_data.Ue.values[0], 
                real_y_locs / bl_data.delta_99.values[0],
                ls='',
                **plot_config
            )

            ax_Uc.plot(
                array(Um) / bl_data.Ue.values[0], 
                real_y_locs / bl_data.delta_99.values[0],
                ls='--',
                color = color,
            )

        ax_Uc.set_xlabel( r"$u_c/u_e$" )
        ax_Uc.set_ylabel( r"$y/\delta_{{99}}$" )
        ax_Uc.set_xlim(0,1.2)
        ax_Uc.set_ylim(0,1.6)

        if schematic:
            im = plt.imread( get_sample_data( schematic  ) )
            newax = fig_Uc.add_axes([0.175, 0.55, 0.3, 0.3], anchor = 'SW', 
                                         zorder=100)
            newax.imshow(im)
            newax.axis('off')

        fig_Uc.savefig(
            "Results/Uc_{0}_{1}.png".format(
                plot_name.replace('.','_').replace('_png','.png'),
                var,
            ), bbox_inches = 'tight'
        )
        fig_std.savefig(
            "Results/std_{0}_{1}.png".format(
                plot_name.replace('.','_').replace('_png','.png'),
                var,
            ), bbox_inches = 'tight'
        )

    for y_loc in coherence_df.near_y_downwind.unique():

        y_df = coherence_df[ coherence_df.near_y_downwind == y_loc ]

        for var in ['u','v']:

            fig_phi,  ax_phi  = plt.subplots(1, 1, figsize = figsize)
            fig_coh,  ax_coh  = plt.subplots(1, 1, figsize = figsize)
            fig_quad, ax_quad = plt.subplots(1, 1, figsize = figsize)

            for case in y_df.case.unique():

                case_df = y_df[ y_df.case == case ]

                bl_data = get_bl_parameters(case_df)
                
                color, marker, cmap = get_color_and_marker( case )

                plot_config = {
                    'marker'          : marker,
                    'markeredgewidth' : markeredgewidth,
                    'markerfacecolor' : markerfacecolor,
                    'markeredgecolor' : color,
                    'markersize'      : markersize,
                    'mew'             : mew,
                    'color'           : color,
                }

                ax_phi.plot(
                    get_Strouhal( case_df['f_'+var], bl_data.delta_99.values[0],
                                 bl_data.Ue.values[0]),
                    case_df['phi_'+var],
                    ls = '',
                    **plot_config
                )
                ax_phi.plot(
                    get_Strouhal( case_df['f_'+var], bl_data.delta_99.values[0],
                                 bl_data.Ue.values[0] ),
                    case_df['f_'+var]*case_df['slope_'+var],
                    '--',
                    lw = 3,
                    color = color,
                )

                ax_coh.plot(
                    case_df['phi_'+var],
                    case_df['gamma_'+var],
                    ls = '',
                    **plot_config
                )

                eta = -0.22
                ax_coh.plot(
                    linspace(0,2*pi,30),
                    exp(eta * linspace(0,2*pi,30)),
                    '--',
                    color = 'k',
                )


            # Configure the phi plot ###########################################
            ax_phi.set_yticks(array(
                [0,1/4.,1/2.,3/4.,1,5./4.,3/2.,7/4.,2]
            )*pi)
            ax_phi.set_yticklabels(
                ['$0$','$\\pi/4$','$\\pi/2$','$3\\pi /4$','$\\pi$',
                 '$5\\pi/4$','$3\\pi/2$','$7\\pi /4$','$2\\pi$'
                ]
            )
            ax_phi.set_xlim(0,St_max)
            ax_phi.set_ylim(0,2*pi)

            ax_phi.set_xlabel(r"$\textrm{{St}}=f\delta/u_e$")
            ax_phi.set_ylabel(
                r"$\phi_{0}$ [rad]"\
                .format(var)
            )
            ax_phi.text(
                0.15,5, 
                r'$y/\delta_{{99}} = {0}$'.format(
                    case_df.near_y_downwind.unique()[0]/delta/1000.,
                ))
            # Save the Phi plot ################################################
            fig_phi.savefig(
                "Results/Phi_{0}_y{1:.1f}_{2}.png".format(
                    plot_name.replace('.','_').replace('_png','.png'),
                    case_df.near_y_downwind.unique()[0]/delta/1000.,
                    var,
                ), bbox_inches = 'tight'
            )
            # ##################################################################

            # Configure the coherence plot #####################################
            ax_coh.set_xticks(array(
                [0,1/2.,1,3/2.,2]
            )*pi)
            ax_coh.set_xticklabels(
                ['$0$', '$\\pi/2$', '$\\pi$', '$3\\pi/2$', '$2\\pi$' ]
            )
            ax_coh.set_ylim(0,1)
            ax_coh.set_xlim(0,2*pi)
            ax_coh.set_xlabel(r"$\phi = \mu_{{x0}}\Delta x$")
            ax_coh.set_ylabel(r"$\gamma_{{{0}}}$".format(var))
            ax_coh.text(
                x = 4,
                y = 0.45,
                ha = 'left',
                s = r'$ e^{{{0:.2f} \phi}}$'.format(eta)
            )
            ax_coh.text(
                pi,0.87,
                r'$y/\delta_{{99}} = {0}$'.format(
                    case_df.near_y_downwind.unique()[0]/delta/1000.,
                ))
            # Save it ##########################################################
            fig_coh.savefig(
                "Results/Coherence_{0}_y{1:.1f}_{2}.png".format(
                    plot_name.replace('.','_').replace('_png','.png'),
                    case_df.near_y_downwind.unique()[0]/delta/1000.,
                    var,
                ), bbox_inches = 'tight'
            )
            plt.cla()


# Constants ####################################################################

nperseg         = 2**6
fs              = 10000

delta           = 9.e-3
Ue              = 20.

St_min          = 0.225
St_max          = 2.4

markeredgewidth = 2
markerfacecolor = 'none'
markersize      = 12
mew             = 4 # Marker edge width

figsize         = (8,7)

# ##############################################################################

