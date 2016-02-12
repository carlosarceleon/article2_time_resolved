#root = "/media/carlos/6E34D2CD34D29783/2015-02_SerrationPIV/TR_Data/"
import seaborn as sns
from matplotlib import rc

St_min = 0.225
St_max = 2.4
delta = 13.9
U     = 20
line_styles = ['--','--','-.',':','--']
markers = [
        #u'o', u'v', u'^', u'<', u'>', u'8', u's', u'p', u'*', u'h', u'H', u'D', u'd'
    'D','s','o','v'
]
markers_full = [
    'o', 'v', '^', '<', '>', '8', 's', 'p', '*', 'h', 'H', 'D', 'd'
]

rc('text',usetex=True)
rc('font',weight='normal')

sns.set_context('paper')
sns.set(font='serif',font_scale=1.5,style='whitegrid')
rc('font',family='serif', serif='cm10')

component_dict = {
    'vx' : "u",
    'vy' : "v",
    'vw' : "w",
}

#cases = [
#    "Slit20R21_a0_p0_U20_z00_tr.h5"  ,
#    #"Sr20R21_a0_p0_U20_z05_tr.h5"    ,
#    #"Slit20R21_a0_p0_U20_z05_tr.h5"  ,
#    "Sr20R21_a0_p0_U20_z10_tr.h5"    ,
#    "Slit20R21_a0_p0_U20_z10_tr.h5"  ,
#    "Sr20R21_a12_p0_U20_z00_tr.h5"   ,
#    "Slit20R21_a12_p0_U20_z00_tr.h5" ,
#    #"Sr20R21_a12_p0_U20_z05_tr.h5"   ,
#    #"Slit20R21_a12_p0_U20_z05_tr.h5" ,
#    "Sr20R21_a12_p0_U20_z10_tr.h5"   ,
#    "Slit20R21_a12_p0_U20_z10_tr.h5" ,
#    "STE_a0_p0_U20_z00_tr.h5"        ,
#    "Sr20R21_a0_p0_U20_z00_tr.h5"    ,
#    "STE_a12_p0_U20_z00_tr.h5"       ,
#]

def get_Strouhal(f,delta,U):
    return f*delta/U

def build_all_cases_popular_lines(
    root='/media/carlos/6E34D2CD34D29783/2015-02_SerrationPIV/TR_Data_NewProcessing'):
    import os
    cases = [c for c in os.listdir(root) if c.endswith('.hdf5')]
    for c in cases:
        print c
        build_popular_line_matrix(
            root = root,
            case = c,
            save_folder='/home/carlos/Documents/PhD/Articles/Conference_AIAASciTech2016/Scripts/time_resolved_scripts/point_data')

def to_db(Pxx):
    from numpy import log10,array
    return array(10*log10(Pxx))

def downsample_db_values(df,vmin,vmax,delta):
    from numpy import arange

    levels = arange(vmin, vmax, delta)

    for i in range(len( levels )):
        if not levels[i]==levels[-1]:
            df.dB.ix[ (df.dB >= levels[i]) & (df.dB < levels[i+1]) ] = \
                    (levels[i]+levels[i+1])/2.

    return df

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

def plot_cross_correlation_locations(
    cases               = [],
    case_names          = [],
    root                = '.',
    x_locs              = [0],
    y_locs              = [0], # Delta normalized
    component           = 'vy',
    plot_name           = 'Correlation_test.png',
    presentation        = True,
    test                = False,
    straight_only_at_TE = True
):

    """ Takes the cases, and plots the crosscorrelation for the 
    requested locations, each location on a new figure
    
    Input:
        case: case name to compare to file names
        root_folder: where to find the pickled point time series
        x_locs: the streamwise locations to plot
        y_locs: the wall-normal locations to plot
        component: the velocity component PSD to plot
        plot_name
    Output:
        Figure
    """
    import matplotlib.pyplot as plt
    import pandas as pd
    from numpy import argmin,array,abs,arctan,sqrt,exp,linspace
    from numpy.random import rand
    from scipy.signal import csd
    import os
    from math import pi
    import matplotlib as mpl

    if presentation:
        rc('font',family='sans-serif', serif='sans-serif')
        mpl.rcParams['text.latex.preamble'] = [
            r'\usepackage{siunitx}'  ,
            r'\sisetup{detect-all}'  ,
            r'\usepackage{sansmath}' ,
            r'\sansmath'               
        ]

    if not len(case_names):
        for c in cases:
            case_names.append(c.replace("_",'-'))

    freq_lower_limit = 300

    def remove_angle_jumps(df):
        from numpy import sign

        df.Phi.loc[df.Phi<0] = \
                df.Phi.loc[df.Phi<0] + pi

        for ix in range(len(df))[:-2]:
            dif = df.Phi.ix[ix+1] - df.Phi.ix[ix]
            if abs(dif) > pi*0.4:
                df.Phi.ix[ix+1] = df.Phi.ix[ix+1] - sign(dif) * pi

        df.Phi.loc[df.Phi<0] = \
                df.Phi.loc[df.Phi<0] + pi

        df.Phi.loc[df.f == df.f.max()] = \
                df.Phi.loc[df.f == df.f.max()] + 2*pi

        return df

    def calculate_Uc(df,delta_x):
        #from scipy.interpolate import interp1d
        from scipy.stats import linregress
        #from numpy import linspace
        #df = pd.DataFrame( data = {
        #    'Phi':Phi,
        #    'f'  :f
        #})
        df = df.sort('f',ascending=True).reset_index(drop=True)

        r_value  = 0
        consider = len(df)
        while r_value**2<0.99:
            df = df.ix[:consider].reset_index(drop=True)
            slope, intercept, r_value, p_value, std_err = linregress(
                df.Phi,
                df.f
            )
            consider -= 1
        
        Uc = 2*pi*slope*delta_x/1000.

        return Uc, df, intercept, slope

    fig_Uc,axes_Uc = plt.subplots(
        len(x_locs),len(y_locs),figsize=(10,10),
        sharex=True,sharey=True
    )
    fig_Phi,axes_Phi = plt.subplots(
        len(x_locs),len(y_locs),figsize=(10,10),
        sharex=True,sharey=True
    )
    fig_Coh,axes_Coh = plt.subplots(
        len(x_locs),len(y_locs),figsize=(10,10),
        sharex=True,sharey=True
    )

    step = 1

    tooth_length = 40.

    for case_name,c_cnt,case_label,marker \
            in zip(cases,range(len(cases)),case_names,
                   markers_full[:len(cases)]):
        if 'a0' in case_name:
            delta = 9.6/1000.
        elif 'a12' in case_name:
            delta = 13.7/1000.
        else:
            delta = 0

        # Build the data frame from pickled data if it's not provided
        print "   Loading {0}".format(case_name)
        case_df = pd.read_hdf(
            os.path.join( root, case_name+"_WallNormalData.hdf5"),
            case_name
        )

        # Normalize the y coordinates to the boundary layer size
        case_df.y = case_df.y*tooth_length/(delta*1000)

        # Get the available coordinates
        df_x_coords = array(sorted(case_df.x.unique(),reverse=False))
        available_x_locs      = []
        available_x_neighbors = []
        for x in x_locs:

            if "STE" in case_name and straight_only_at_TE:
                x = min(x_locs)

            x_av, dx = find_nearest(x, df_x_coords)
            
            neighbor_index = argmin(abs(df_x_coords-x_av))+step

            if neighbor_index < len(df_x_coords):
            
                available_x_neighbors.append(
                    df_x_coords[neighbor_index]
                )
                available_x_locs.append(x_av)

        if len(available_x_locs) and len(available_x_neighbors):
            for x_l,x_n,xi in zip(
                available_x_locs, available_x_neighbors,
                range(len(available_x_locs))
            ):

                df_y_coords          = \
                        case_df[case_df.x==x_l].y.unique()
                df_y_coords_neighbor = \
                        case_df[case_df.x==x_n].y.unique()

                for y_l,yi in zip(y_locs,
                                  range(len(y_locs))):

                    y_av, yd = find_nearest(y_l , df_y_coords)
                    y_n,  yd = find_nearest(y_l , df_y_coords_neighbor)

                    plt_idx = len(y_locs)-yi-1

                    time_series = case_df[
                        (case_df.x == x_l) &\
                        (case_df.y == y_av)
                    ].sort('ti').reset_index(drop=True)

                    time_series_neighbor = case_df[
                        (case_df.x == x_n) &\
                        (case_df.y == y_n)
                    ].sort('ti').reset_index(drop=True)
                    
                    non_null_time_series =  time_series[
                        time_series.vx.notnull()
                    ]
                    non_null_time_series_neighbor =  \
                            time_series_neighbor[
                        time_series_neighbor.vx.notnull()
                    ]

                    text = axes_Phi[plt_idx][xi].text( 
                        x = 0.90,
                        y = 0.10,
                        s = "$x/2h = {0:.1f}$, $y/\\lambda = {1:.1f}$"\
                        .format(x_l,y_l),
                        ha = 'right',
                        transform = axes_Phi[plt_idx][xi].transAxes,
                        zorder = 10
                    )
                    text.set_bbox(dict(color='white', alpha=0.5))
                    text = axes_Coh[plt_idx][xi].text( 
                        x = 0.90,
                        y = 0.10,
                        s = "$x/2h = {0:.1f}$, $y/\\lambda = {1:.1f}$"\
                        .format(x_l,y_l),
                        ha = 'right',
                        transform = axes_Coh[plt_idx][xi].transAxes,
                        zorder = 10
                    )
                    text.set_bbox(dict(color='white', alpha=0.5))
                    text = axes_Uc[plt_idx][xi].text( 
                        x = 0.10,
                        y = 0.10,
                        s = "$x/2h = {0:.1f}$, $y/\\lambda = {1:.1f}$"\
                        .format(x_l,y_l),
                        transform = axes_Uc[plt_idx][xi].transAxes,
                        zorder = 10
                    )
                    text.set_bbox(dict(color='white', alpha=0.5))

                    if not non_null_time_series.empty\
                       and not non_null_time_series_neighbor.empty\
                       and len(non_null_time_series_neighbor) == \
                       len(non_null_time_series):
                        
                        max_lag = 10000
                        s1 = non_null_time_series[component]\
                                .values[0:max_lag] \
                                - non_null_time_series[component]\
                                .values[0:max_lag].mean()
                        s2 = non_null_time_series_neighbor[component]\
                                .values[0:max_lag] \
                                - non_null_time_series_neighbor[component]\
                                .values[0:max_lag].mean()

                        if test:
                            s1 = rand(max_lag)
                            s2 = rand(max_lag)

                        f,Pxy = csd(
                            s2,s1,
                            nperseg = 2**6,
                            fs      = 10000,
                        )

                        f,Pxx = csd(
                            s1,s1,
                            nperseg = 2**6,
                            fs      = 10000,
                        )

                        f,Pyy = csd(
                            s2,s2,
                            nperseg = 2**6,
                            fs      = 10000,
                        )

                        gamma_squared = \
                                abs(Pxy)**2 / ( Pxx * Pyy )

                        gamma = sqrt(gamma_squared)

                        Phi = arctan( Pxy.imag / Pxy.real )

                        df = pd.DataFrame( data = {
                            'Phi':Phi,
                            'f'  :f,
                            'gamma': gamma
                        })

                        df = df[df.f >= freq_lower_limit].reset_index(
                            drop = True
                        )

                        df = remove_angle_jumps(df)
                        df = remove_angle_jumps(df)

                        line = axes_Phi[plt_idx][xi].plot(
                            get_Strouhal(df.f,delta,U),
                            df.Phi,
                            alpha = 0.3,
                        )


                        eta = 0.22
                        axes_Coh[plt_idx][xi].plot(
                            linspace(0,2*pi,30),
                            exp(-eta * linspace(0,2*pi,30)),
                            '--',
                            color = 'k',
                        )

                        axes_Coh[plt_idx][xi].scatter(
                            df.Phi,
                            df.gamma,
                            color  = line[0].get_color(),
                            alpha  = 0.3,
                            marker = marker
                        )

                        Uc,df,intercept,slope = calculate_Uc(
                            df,
                            delta_x = abs(x_n - x_l) * tooth_length
                        )

                        df.Strouhal = get_Strouhal(df.f,delta,U)

                        axes_Coh[plt_idx][xi].scatter(
                            df.Phi,
                            df.gamma,
                            color = line[0].get_color(),
                            label = case_label,
                            marker = marker
                        )

                        axes_Phi[plt_idx][xi].plot(
                            df.Strouhal,
                            df.Phi,
                            color = line[0].get_color(),
                            label = case_label
                        )

                        if df.f.max()>1000:

                            axes_Phi[plt_idx][xi].plot(
                                df.Strouhal,
                                df.f*slope**(-1),
                                '--',
                                color = line[0].get_color(),
                            )

                            bar_width = 1.
                            axes_Uc[plt_idx][xi].bar(
                                left   = c_cnt+bar_width/2.5,
                                width  = bar_width*0.8,
                                color  = line[0].get_color(),
                                height = Uc/20.,
                                label  = case_label
                            )

    
    for axi in axes_Phi:
        for ax in axi:
            ax.set_yticks(array(
                [0,1/4.,1/2.,3/4.,1,5./4.,3/2.,7/4.,2]
            )*pi)
            ax.set_yticklabels(
                ['$0$','$\\pi/4$','$\\pi/2$','$3\\pi /4$','$\\pi$',
                 '$5\\pi/4$','$3\\pi/2$','$7\\pi /4$','$2\\pi$'
                ]
            )
            ax.set_xlim(St_min,St_max)
            ax.set_ylim(0,2*pi)
            ax.set_xlabel("")
            ax.set_ylabel("")

    for axi in axes_Coh:
        for ax in axi:
            ax.set_xticks(array(
                [0,1/2.,1,3/2.,2]
            )*pi)
            ax.set_xticklabels(
                ['$0$','$\\pi/2$','$\\pi$',
                 '$3\\pi/2$','$2\\pi$'
                ]
            )
            ax.set_ylim(0,1)
            ax.set_xlim(0,2*pi)
            ax.set_xlabel("")
            ax.set_ylabel("")

    for axi in axes_Uc:
        for ax in axi:
            #ax.set_xscale('log')
            ax.set_xticks(range(len(cases)))
            ax.set_xticklabels(['']*len(cases))
            ax.set_ylim(0.2,1.2)
            ax.set_xlabel("")
            ax.set_ylabel("")


    axes_Phi[len(x_locs)-1][0].\
            set_xlabel("$\\textrm{{St}}_\\delta$")
    axes_Phi[len(x_locs)-1][0].\
            set_ylabel(
                "$\\phi_{{x,x+\\Delta x}},\, {0}$ [rad]"\
                .format(component_dict[component],x_n-x_l)
            )
    axes_Coh[len(x_locs)-1][0].\
            set_xlabel("$\\phi = \mu_{{x0}}\\Delta x$")
    axes_Coh[len(x_locs)-1][0].\
            set_ylabel(
                "$\\gamma$"
            )
    axes_Uc[len(x_locs)-1][0].\
            set_xlabel("$\\textrm{{St}}_\\delta$")
    axes_Uc[len(x_locs)-1][0].\
            set_ylabel(
                "$U_{{c}}/U_{{\infty}}$"\
                .format(component_dict[component])
            )

    axes_Phi[0][0].legend(
        bbox_to_anchor = (0., 1.02, len(x_locs), .102),
        loc            = 3,
        ncol           = 2,
        mode           = "expand",
        borderaxespad  = 0.
    )
    axes_Coh[0][0].legend(
        bbox_to_anchor = (0., 1.02, len(x_locs), .102),
        loc            = 3,
        ncol           = 2,
        mode           = "expand",
        borderaxespad  = 0.
    )
    axes_Uc[0][0].legend(
        bbox_to_anchor = (0., 1.02, len(x_locs), .102),
        loc            = 3,
        ncol           = 2,
        mode           = "expand",
        borderaxespad  = 0.
    )

    axes_Coh[0][0].annotate(
        "$\\textrm{exp}\\left(-\\eta\\phi\\right)$", 
        xy=(pi/2., exp(-eta*pi/2.)), xycoords='data',
        xytext=(pi/2.+pi/2., exp(-eta*pi/2.)+0.1), 
        textcoords='data',
        size=15,
        # bbox=dict(boxstyle="round", fc="0.8"),
        arrowprops=dict(
            arrowstyle='simple',
            fc="k", ec="w",
            #patchB=el,
            connectionstyle="arc3,rad=0.3",
        ),
    )

    fig_Phi.savefig(
        plot_name.replace('.png','_PhaseSpectra.png'), 
        bbox_inches='tight'
    )
    fig_Coh.savefig(
        plot_name.replace('.png','_Coherence.png'), 
        bbox_inches='tight'
    )
    fig_Uc.savefig(
        plot_name.replace('.png','_ConvectionVelocity.png'), 
        bbox_inches='tight'
    )
    return 0

def plot_select_frequency_locations(
    cases               = [],
    case_names          = [],
    root                = '.',
    x_locs              = [0],
    y_locs              = [0], # Delta normalized
    component           = 'vy',
    plot_name           = 'PSD_test.png',
    scale_by_freq       = False,
    presentation        = True,
    straight_only_at_TE = True,
):

    """ Takes the cases, and plots the PSD for the requested locations, 
    each location on a new figure
    
    Input:
        case: case name to compare to file names
        root_folder: where to find the pickled point time series
        x_locs: the streamwise locations to plot
        y_locs: the wall-normal locations to plot
        component: the velocity component PSD to plot
        plot_name
        scale_by_freq
    Output:
        Figure
    """
    import matplotlib.pyplot as plt
    import pandas as pd
    from numpy import linspace,log10,arange
    import os
    from scipy.signal import welch
    import matplotlib as mpl

    if presentation:
        rc('font',family='sans-serif', serif='sans-serif')
        mpl.rcParams['text.latex.preamble'] = [
            r'\usepackage{siunitx}'  ,
            r'\sisetup{detect-all}'  ,
            r'\usepackage{sansmath}' ,
            r'\sansmath'               
        ]

    if not len(case_names):
        for c in cases:
            case_names.append(c.replace("_",'-'))

    fig_vel,axes_vel = plt.subplots(
        len(x_locs),len(y_locs),figsize=(10,10),
        sharex=True,sharey=True
    )
    fig_freq,axes_freq = plt.subplots(
        len(x_locs),len(y_locs),figsize=(10,10),
        sharex=True,sharey=True
    )

    tooth_length = 40.
    for case_name,case_label in zip(cases,case_names):
        if 'a0' in case_name:
            delta = 9.6/1000.
        elif 'a12' in case_name:
            delta = 13.7/1000.
        else:
            delta = 0

        slope_x = linspace(0.7,1.5,100)
        if component == 'vx':
            vmax = 0
            vmin = -25
            yticks = arange(vmin,vmax,5)
            slope_y = 10*log10(slope_x**(-5/3.))-7
            f_loc = -6

        elif component == 'vy':
            vmax = 0
            vmin = -25
            yticks = arange(vmin,vmax,5)
            slope_y = 10*log10(slope_x**(-5/3.))-7 - 4
            f_loc = -6 - 4

        # Build the data frame from pickled data if it's not provided
        print "   Loading {0}".format(case_name)
        case_df = pd.read_hdf(
            os.path.join( root, case_name+"_WallNormalData.hdf5"),
            case_name
        )

        # Normalize the y coordinates to the boundary layer size
        case_df.y = case_df.y*tooth_length/(delta*1000)

        # Get the available coordinates
        df_x_coords = sorted(case_df.x.unique())
        available_x_locs = []
        for x in x_locs:

            if "STE" in case_name and straight_only_at_TE:
                x = min(x_locs)

            x_av, dx = find_nearest(x, df_x_coords)
            available_x_locs.append(x_av)

        for x_l,xi in zip(available_x_locs,
                          range(len(available_x_locs))):
            # Create the dataframe that will hold the wall
            # normal PSD data

            df_y_coords = case_df[case_df.x==x_l].y.unique()

            for y_l,yi in zip(y_locs,
                              range(len(y_locs))):

                y_av, yd = find_nearest(y_l,df_y_coords)

                time_series = case_df[
                    (case_df.x == x_l) &\
                    (case_df.y == y_av)
                ].sort('ti')
                
                non_null_time_series =  time_series[
                    time_series.vx.notnull()
                ]
                
                text_freq = axes_freq[len(y_locs)-yi-1][xi].text( 
                    x = 0.10,
                    y = 0.10,
                    ha = 'left',
                    s="$x/2h = {0:.1f}$, $y/\\lambda = {1:.1f}$".format(x_l,y_l),
                    transform = axes_freq[len(y_locs)-yi-1][xi]\
                    .transAxes,
                    zorder=10
                )
                text_freq.set_bbox(dict(color='white', alpha=0.5))

                text_vel = axes_vel[len(y_locs)-yi-1][xi].text( 
                    x = 0.90,
                    y = 0.90,
                    s=\
                    "$x/2h = {0:.1f}$, $y/\\lambda = {1:.1f}$"\
                    .format(x_l,y_l),
                    transform = axes_vel[len(y_locs)-yi-1][xi]\
                    .transAxes,
                    zorder=10
                )
                text_vel.set_bbox(dict(color='white', alpha=0.5))

                if not non_null_time_series.empty:

                    if case_name == 'Slit20R21_phi0_alpha12_U20_loc10':
                        Fs = 5000
                    else:
                        Fs = 10000
                    #fig_tmp,ax_tmp = plt.subplots(1,1)
                    #Pxx,freq = ax_tmp.psd(
                    #    non_null_time_series[component].values,
                    #    NFFT          = 2**6,
                    #    Fs            = Fs,
                    #    scale_by_freq = scale_by_freq,
                    #)
                    #plt.close(fig_tmp)
                    freq, Pxx = welch(
                        non_null_time_series[component].values,
                        nperseg          = 2**7,
                        fs            = Fs,
                        scaling = 'spectrum',
                    )

                    axes_freq[len(y_locs)-yi-1][xi].plot(
                        get_Strouhal(freq,delta,U),
                        10*log10(Pxx),
                        label = case_label
                    )
                    non_null_time_series[non_null_time_series.ti<1000]\
                            .plot(
                                x     = 't_real',
                                y     = component,
                                ax    = axes_vel[len(y_locs)-yi-1][xi],
                                label = case_label,
                            )
                    return_df = non_null_time_series

        for axi in axes_vel:
            for ax in axi:
                ax.set_xlim(0,1000/10000.)
                ax.set_ylim(-1,25)
                ax.set_ylabel(component_dict[component]+" [m/s]")
                ax.set_xlabel("Time [s]")
                ax.set_xlabel("")
                ax.set_ylabel("")

        for axi in axes_freq:
            for ax in axi:
                ax.set_xscale('log')
                #ax.set_xticks(Strouhal_range)
                #ax.set_xticklabels(
                #    ['0.15', '0.5', '1', '1.5', '2', '2.5']
                #)
                ax.set_yticks(yticks)
                #ax.grid(False)
                ax.set_xlim(St_min,St_max)
                ax.set_ylim(vmin, vmax)
                ax.plot(slope_x,slope_y,color='k',lw=1)
                ax.set_xlabel("")
                ax.set_ylabel("")

        axes_freq[len(x_locs)-1][0].\
                set_ylabel(
                    "$10\\log_{{10}}\\left(\\Phi_{{\\textrm{{TKE}},{0}}}\\right)$ [dB]".format(
                        component_dict[component])
                )
        axes_freq[len(x_locs)-1][0].\
                set_xlabel("$\\textrm{{St}}_\\delta$")

        axes_freq[0][0].text(1.0,f_loc,"$\\textrm{St}^{-5/3}_\\delta$",fontsize=15)

        axes_vel[0][0].legend(
            bbox_to_anchor = (0., 1.02, len(x_locs), .102),
            loc            = 3,
            ncol           = 2,
            mode           = "expand",
            borderaxespad  = 0.
        )
        axes_freq[0][0].legend(
            bbox_to_anchor = (0., 1.02, len(x_locs), .102),
            loc            = 3,
            ncol           = 2,
            mode           = "expand",
            borderaxespad  = 0.
        )

    #axes[0].set_ylabel("$y/\\delta_{95}$")
    #axes[1].set_xlabel("$f$ [kHz]")

    fig_vel.savefig(
        plot_name.replace('.png','_Velocity.png'), bbox_inches='tight'
    )
    fig_freq.savefig(
        plot_name.replace('.png','_PSD.png'), bbox_inches='tight'
    )
    return return_df

def make_wall_normal_frequency_map(case_name,root,
                                   plot_name='FrequencyMap.png',
                                   df=None,
                                   presentation=False):
    """ Takes the case, and looks for the available data points in the
    root directory. Builds a vertical (same x) frequency heat map
    from the time series found in there

    Input:
        case: case name to compare to file names
        root_folder: where to find the pickled point time series
        plot_name
        df: optional, the data frame with the time series
    Output:
        Figure
    """
    from numpy import meshgrid
    import matplotlib.pyplot as plt
    import pandas as pd
    import matplotlib as mpl
    from matplotlib import rc
    from os.path import split
    from matplotlib import mlab

    #import seaborn as sns

    print plot_name
    if presentation:
        rc('font',family='sans-serif', serif='sans-serif')
        mpl.rcParams['text.latex.preamble'] = [
            r'\usepackage{siunitx}',   # i need upright \micro symbols, but you need...
            r'\sisetup{detect-all}',   # ...this to force siunitx to actually use your fonts
            r'\usepackage{sansmath}',  # load up the sansmath so that math -> helvet
            r'\sansmath'               # <- tricky! -- gotta actually tell tex to use!
        ]

    scale_by_freq = False
    if scale_by_freq:
        vmin = -70
        vmax = -25
    else:
        vmin = -70+35
        vmax = -25+35

    if 'a0' in case_name:
        delta = 9.6/1000.
    elif 'a12' in case_name:
        delta = 13.7/1000.
    else:
        delta = 0

    # Build the data frame from pickled data if it's not provided
    if df is None:
        case_df = pickled_coordinate_files_to_DF(
            case_name = case_name,
            root      = root
        )
    else:
        case_df = df

    # Get the available coordinates
    x_coords = sorted(case_df.x.unique())
    for x in x_coords:
        y_coords = case_df[case_df.x==x].y.unique()

    if len(x_coords)==1:
        extra_frame = 1
    else:
        extra_frame = 0
    fig,ax = plt.subplots(1,len(x_coords)+extra_frame,figsize=(20,5),
                         sharex=True,sharey=True)
    for x,xi in zip(x_coords,range(len(x_coords))):
        # Create the dataframe that will hold the wall
        # normal PSD data
        db_df = pd.DataFrame(columns=['y','freq','dB'])

        for y in y_coords:
            time_series = case_df[
                (case_df.x == x) &\
                (case_df.y == y)
            ].sort('ti')

            fig2,ax2 = plt.subplots(1,1)
            if case_name == 'Slit20R21_phi0_alpha12_U20_loc10':
                Fs = 5000
            else:
                Fs = 10000
            vel_component = time_series.vx.values.astype('float')
            Pxx,freq = ax2.psd(
                x             = vel_component**2,
                NFFT          = 2**6,
                Fs            = Fs,
                window        = mlab.window_hanning,
                scale_by_freq = scale_by_freq)
            plt.close(fig2)
            

            db_df = db_df.append(
                pd.DataFrame( data = {
                    'y'    : [y*40]*len(freq),
                    'freq' : freq,
                    'dB'   : to_db(Pxx)
                })
            ).sort(['y','freq']).fillna(0)

        # Ignore all frequencies above or equal 200 Hz
        db_df = db_df[db_df.freq >= 200]
        db_df = db_df[db_df.dB > vmin]
        db_df.y = db_df.y/(delta*1000)
        db_df.freq = db_df.freq/1000.
        #db_df = downsample_db_values(db_df,vmin,vmax,2.5)
        X,Y = meshgrid(db_df.freq.unique(),db_df.y.unique())
        try:
            Z   = db_df.dB.reshape(X.shape)
        except ValueError:
            print db_df
            return 0

        # Turn the dataframe into a pivot table, for the heatmap
        db_df_pivot = db_df.pivot(index='y',columns='freq',values='dB')

        # Turn it the right side up (for the y axis)
        db_df_pivot = db_df_pivot.sort_index(axis=0,ascending=False)

        fmap = ax[xi].pcolor(
            X,Y,Z,
            vmin = vmin,
            vmax = vmax,
            cmap = 'RdBu_r',
            zorder = 5
            #interpolation = 'nearest'
        )
        ax[xi].set_xlim(0.250, 5)
        ax[xi].set_ylim(0, 2)
        text = ax[xi].text( x=0.62,y=0.87,s="$x/2h = {0:.1f}$"\
                        .format(x),
                        fontsize=25,
                        transform = ax[xi].transAxes,
                        zorder=10
            
        )
        text.set_bbox(dict(color='white', alpha=0.5))

        #ax[xi].contour(X,Y,Z)
        ax[xi].set_xscale('log')
        ax[xi].set_xticks([1,2,3,4,5])
        ax[xi].set_xticklabels(['1','2','3','4','5'])
        ax[xi].grid(False)
    ax[0].set_ylabel("$y/\\delta_{95}$")
    ax[1].set_xlabel("$f$ [kHz]")
    ax[1].set_title(case_name)
    fig.subplots_adjust(right=0.88)
    cbar_ax = fig.add_axes([0.90, 0.15, 0.02, 0.7])
    clb = fig.colorbar(fmap, cax=cbar_ax)
    clb.set_label("PSD [dB]",labelpad=20)
    if presentation:
        plt.savefig(
            split(plot_name)[0]+"/Presentation_"+split(plot_name)[1],
                                                    bbox_inches='tight')
    else:
        plt.savefig(plot_name, bbox_inches='tight')
        

def test_NFFT(s):
    from numpy import arange,power
    from pylab import psd,show
    import matplotlib.pyplot as plt 
    import seaborn as sns
    sns.__version__

    NFFTs = power(2,arange(6,10))

    for NFFT in NFFTs[::-1]:
        psd(s,NFFT,scale_by_freq=True,Fs=10000)

    plt.xlim(200,5000)
    plt.ylim(-45,-35)
    plt.semilogx()
    show()


def test_psd():

    from numpy import linspace,sin
    from pylab import psd,show
    x = linspace(0,100,10000)
    y = sin(x*100)/100.

    psd(y,NFFT=256,scale_by_freq=True)
    show()

def all_pickled_coordinates_to_DFs(root,overwrite=False):
    from re import findall
    from numpy import array,unique
    from os import listdir,path

    available_case_point_files = [f for f in listdir(root)\
                                  if f.endswith('.p')]

    cases = []
    for available in available_case_point_files:
        cases.append(
            findall("[A-Za-z0-9_]+_px",available)[0].replace("_px","")
        )
    cases = unique(array(cases))

    print "Found the following cases, which will turn into data frames"
    for c in cases:
        print "\t{0}".format(c)

    if not overwrite:
        cases_to_run = []
        for c in cases:
            if not path.isfile(c+".p"):
                cases_to_run.append(c)
        cases = cases_to_run
        print "Will run the following non-processed cases"
        for c in cases:
            print "\t{0}".format(c)

    DFs = []
    for case in cases:
        DFs.append(pickled_coordinate_files_to_DF(case,root))

    return DFs

def pickled_coordinate_files_to_DF(case_name,root):
    """ Reads the files in the root that contain the given case name
    and compiles a data frame with all the available coordinate
    timeseries

    Input:
        case_name: case name
        root:      where to find the pickled case coordinate points

    Output:
        pandas data frame
    """

    import os
    import pandas as pd
    from progressbar import ProgressBar,Percentage,Bar
    from progressbar import ETA,SimpleProgress

    available_case_point_files = [f for f in os.listdir(root)\
                                  if case_name in f]

    print "Loading case coordinates into data frame..."
    print case_name
    progress = ProgressBar(
        widgets=[
            Bar(),' ',
            Percentage(),' ',
            ETA(), ' (file ',
            SimpleProgress(),')'], 
        maxval=len(available_case_point_files)
    ).start()

    case_time_series_df = pd.read_pickle(
        os.path.join(root,available_case_point_files[0])
    )

    cnt_files = 1
    for coord_file in available_case_point_files[1:]:
        case_time_series_df = pd.concat(
            [case_time_series_df,
             pd.read_pickle(os.path.join(root,coord_file))]
        )
        cnt_files+=1
        progress.update(cnt_files)
    progress.finish()

    case_time_series_df['case'] = case_name

    return case_time_series_df

def load_time_series(device="Sr20R21",alpha='0',phi='0',z='10', p=(-0.1,0.1)):
    from time_data_functions import read_hdf5_time_series
    import os

    #case = "{0}_a{1}_p{2}_U20_z{3}_tr_planar".format(
    #        device,
    #        alpha,
    #        phi,
    #        z
    #)
    case = "{0}_a{1}_p{2}_U20_z{3}_tr_NewProcessing".format(
            device,
            alpha,
            phi,
            z
    )

    hdf5_file = os.path.join(
        root,
        'TimeData_NewProcessing.hdf5'
        #'planar.hdf5'
    )
    
    return read_hdf5_time_series(
            hdf5_file,
            case,
            loc = p
        ).interpolate()

def check_PSD(device="Sr20R21",alpha='0',phi='0',z='10',
              p=(-0.1,0.1),component='vy',plot_name='test_psd.png',
             sampling=10000,NFFT=256):
    from time_series_functions import get_spectral_power_density#,\
            #butter_lowpass_filter
    from matplotlib import pyplot as plt
    import seaborn as sns
    sns.__version__

    time_series = load_time_series(device=device,alpha=alpha,phi=phi,z=z,p=p)

    #time_series.vx = butter_lowpass_filter(time_series.vx,cutoff=2000,fs=10000)
    #time_series.vy = butter_lowpass_filter(time_series.vy,cutoff=2000,fs=10000)
    #time_series.vz = butter_lowpass_filter(time_series.vz,cutoff=2000,fs=10000)
    
    Pxx,freqs = get_spectral_power_density(time_series.vx,
                                           NFFT=NFFT,sampling=sampling)
    Pyy,freqs = get_spectral_power_density(time_series.vy,
                                           NFFT=NFFT,sampling=sampling)
    #Pzz,freqs = get_spectral_power_density(time_series.vz,
    #                                       NFFT=NFFT,sampling=sampling)

    if plot_name:
        fig = plt.figure()

        plt.plot(freqs, Pxx, label='$u$', alpha=0.6)
        plt.plot(freqs, Pyy, label='$v$', alpha=0.6)
        #plt.plot(freqs, Pzz, label='$z$', alpha=0.6)

        plt.ylabel("Power spectral density [p$^2/$Hz]")
        plt.xlabel("St")
        plt.yscale('log')
        plt.xscale('log')
        plt.xlim(xmin=50,xmax=5000)
        #plt.ylim(ymin=10e-4)
        plt.legend(loc='upper right')
        plt.savefig(plot_name)
        fig.clear()
    return [Pxx,Pyy,freqs]


def check_autocorrelation(device="Sr20R21",alpha='0',phi='0',z='10',
                          p=(-0.1,0.1),component='vx',plot_name='test.png'):
    from time_series_functions import get_autocorrelation,butter_lowpass_filter
    from matplotlib import pyplot as plt
    import seaborn as sns
    from numpy import arange
    sns.__version__

    time_series = load_time_series(device=device,alpha=alpha,phi=phi,z=z,p=p)

    time_series.vx = butter_lowpass_filter(time_series.vx,cutoff=2000,fs=10000)
    time_series.vy = butter_lowpass_filter(time_series.vy,cutoff=2000,fs=10000)
    time_series.vz = butter_lowpass_filter(time_series.vz,cutoff=2000,fs=10000)

    autocorr_u = get_autocorrelation(
        time_series.vx
    )
    autocorr_v = get_autocorrelation(
        time_series.vy
    )
    autocorr_w = get_autocorrelation(
        time_series.vz
    )

    if plot_name:
        fig = plt.figure()
        plt.plot(arange(len(autocorr_u)),
                 autocorr_u/autocorr_u.max(),label='$u$',alpha=0.6)
        plt.plot(arange(len(autocorr_u)),
                 autocorr_v/autocorr_v.max(),label='$v$',alpha=0.6)
        plt.plot(arange(len(autocorr_u)),
                 autocorr_w/autocorr_w.max(),label='$w$',alpha=0.6)
        plt.xlabel("Lag [time steps (1:1/10000 s)]")
        plt.ylabel("Autocorrelation")
        plt.xlim(0,20) 
        plt.legend(loc='upper right')
        plt.savefig(plot_name)
        fig.clear()

    return autocorr_u,autocorr_v

def plot_time_series(device="Sr20R21",alpha='0',phi='0',z='10',
                     p=(-0.1,0.1),component='vx',plot_name='test.png'):
    from matplotlib import pyplot as plt
    import seaborn as sns
    from time_series_functions import butter_lowpass_filter
    sns.__version__

    time_series = load_time_series(device=device,alpha=alpha,
                                   phi=phi,z=z,p=p)

    
    time_series_low_passed = butter_lowpass_filter(time_series.vx,
                                                   cutoff=2000,fs=10000)

    fig = plt.figure()
    plt.plot(time_series.t,time_series.vx,label='$u$',alpha=0.6)
    plt.plot(time_series.t,time_series_low_passed,
             label='Low passed $u$',alpha=1.0,lw=3,color='k')
    #plt.plot(time_series.t,time_series.vy,label='$v$',alpha=0.6)
    #plt.plot(time_series.t,time_series.vz,label='$w$',alpha=0.6)
    plt.xlabel("t [s]")
    plt.ylabel("Velocity [m/s]")
    plt.xlim(0,500/10000.)
    plt.ylim(-5,22)
    plt.legend(loc='lower right')
    plt.savefig(plot_name)
    fig.clear()

def popular_points():
    from numpy import linspace
    px = linspace(0,1,5)
    # 2h = 4.0cm; the BL is about 1.5 cm, so 0.375*2h, plus ignore
    # the first 0.4 cm, so it starts at 0.1
    py = linspace(0.1,0.375,6)
    return px,py

def popular_lines():
    px = [0,0.5,1]

    return px

def build_popular_point_matrix(root,case,save_folder=0):
    """ Get a 5x5 matrix of the most important (used) points 
    for processing, and save them as pickles of the time series

    Input
        case:           the HDF5 case
        save_folder:    folder where to save the pickle
    """

    if not save_folder:
        save_folder = root

    px,py = popular_points()

    for x in px:
        for y in py:
            make_pickled_time_series(root     = root,
                                  case        = case,
                                  x           = x,
                                  y           = y,
                                  save_folder = save_folder
                                 )

def build_popular_line_matrix(root,case,save_folder=0):
    """ Get a matrix of the most important (used) points 
    for processing, and save them as pickles of the time series

    Input
        case:           the HDF5 case
        save_folder:    folder where to save the pickle
    """
    from numpy import linspace

    if not save_folder:
        save_folder = root

    px = popular_lines()

    for x in px:
        # 2h = 4.0cm; the BL is about 1.5 cm, so 0.375*2h, plus ignore
        # the first 0.4 cm, so it starts at 0.1
        for y in linspace(0.0,1.5,100):
            make_pickled_time_series(root        = root,
                                     case        = case,
                                     x           = x,
                                     y           = y,
                                     save_folder = save_folder,
                                     overwrite   = False
                                 )

def make_pickled_time_series(root,case,x,y,save_folder,overwrite=False):
    import time_data_functions as tdf
    from os.path import join,isfile
    if not overwrite and not isfile(
        join(save_folder,"{0}_px{1}_py{2}.p".format(
            case.replace('.hdf5',''),
            x,
            y
        ))):
        df = tdf.read_hdf5_time_series(
            join(root,case),
            case.replace('.hdf5',''),
            loc=(-y,x)
        )
    else:
        return 0
    if df is not None:
        df['x'] = x
        df['y'] = y
        df.to_pickle(join(save_folder,"{0}_px{1}_py{2}.p".format(
            case.replace('.hdf5',''),
            x,
            y
        )))
        return 0

def check_quick_point_existence(requested_data_points,root,data_target):
    """ This function checks if the points requested exist amongst
    the quick points popular points matrix, otherwise go and make them

    Input: a file array of the available points for this device and
        conditions

    Returns: the files that have been approved for processing

    """
    from numpy import argmin
    from re import findall

    px,py = popular_points()
    data_points  = []
    non_existant = []
    for all_data in requested_data_points:
        x = float(findall('px[0-9]+.[0-9]+',all_data)[0].replace('px',''))
        y = float(findall('py[0-9]+.[0-9]+',all_data)[0].replace('py',''))
        dy = abs(y-py)
        dx = abs(x-px)
        if dx[argmin(dx)]<1e-5 and dy[argmin(dy)]<1e-5:
            data_points.append(all_data)
        else:
            non_existant.append(all_data)

    if len(non_existant):
        for n_existant in non_existant:
            case = findall("^[A-Za-z0-9.-]_tr",n_existant)[0]+".hdf5"
            make_pickled_time_series(root,case,x,y,data_target)
            data_points.append(n_existant)

    return data_points


#p = (-0.05,0.3)
#check_PSD(p=p,NFFT=512)
#acorr = check_autocorrelation(plot_name='test_autocorrelation_new.png')
#plot_time_series(plot_name='test_timeseries_new.png')
#df = plot_surface_at_t(p=p,plot_name='test_surface.png')
#for i in range(100):
#    df = plot_surface_at_t(t=i,plot_name='test_{0:03d}.png'.format(i))

