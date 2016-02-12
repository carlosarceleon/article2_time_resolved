def butter_lowpass(cutoff, fs, order=5):
    from scipy.signal import butter
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    return b, a

def butter_lowpass_filter(data, cutoff, fs, order=5):
    from scipy.signal import lfilter
    from numpy import array
    b, a = butter_lowpass(cutoff, fs, order=order)
    y = lfilter(b, a, data)
    y = [data[0]]+list(y[1:])
    return array(y)

def get_spectral_power_density(s1,sampling=5000,NFFT=256):
    """ Gets the power spectral density of the signal s1

    Input:
        s1: time series
        sampling: the sampling of the signal [def: 5000]
        NFFT: The number of data points used in each block for the FFT

    Output:
        Pxx: the power spectral density
        freqs: the frequencies
    """
    from matplotlib.mlab import psd

    return psd(s1,Fs=sampling,NFFT=NFFT,scale_by_freq=True)


def get_coherence(s1,s2,sampling=5000):
    """ Gets the coherence between signals s1 and s2

    Input:
        s1,s2: two time series
        sampling: the sampling of the signal [def: 5000]

    Output:
        Cxy
    """
    from scipy.signal import coherence

    return coherence(s1,s2,sampling=sampling)

def get_autocorrelation(s1,mode='full'):
    """ Gets the autocorrelation of signal s1

    Input:
        s1: a time series
        mode: see the arguments in scipy.signa.correlate

    Output:
        Ruu
    """
    from scipy.signal import correlate

    return correlate(s1-s1.mean(),s1-s1.mean(),mode=mode)[len(s1)+1:]


def get_crosscorrelation(s1,s2,mode='full'):
    """ Gets the cross-correlation between signals s1 and s2

    Input:
        s1,s2: two time series
        mode: see the arguments in scipy.signa.correlation

    Output:
        Ruu
    """
    from scipy.signal import correlate

    return correlate(s1,s2,mode=mode)

def get_turbulence_length_scale(s1,U=20,mode='full'):
    """ Gets the turbulence length scale of signal s1.
    It assumes that the turbulence is frozen, and thus
                L = U_c * \\tau_0
    where U_c ~ 0.6 U, and \\tau_0 is the value of the
    correlation lag, for which
                R_uu(x,\\tau_0) / R_uu(x,0) = 0.5

    Inputs:
        s1: a time series
        U: the freestream velocity
        mode: see the arguments in scipy.signa.correlate
    """
    from scipy.signal import correlate

    convection_U = 0.6*U

    R_uu = correlate(s1,s1,mode=mode)
    for r,tau in zip(R_uu,range(len(R_uu))):
        if r/R_uu[0] < 0.5:
            break

    return tau*convection_U

def get_coherence_phase(s1,s2,sampling=5000):
    """ Gets the phase angle of the coherence between two signals

    Input:
        s1,s2: two time series
        sampling: the sampling of the signal [def: 5000]

    Output:
        Psi: phase in radians
        f: frequencies
    """
    from matplotlib.mlab import csd
    from numpy import angle

    Pxy, freqs = csd( s1,s2, NFFT = 256, Fs = sampling )
    return angle( Pxy ), freqs

