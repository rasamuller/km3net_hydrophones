# -*- coding: utf-8 -*-
"""
Created on Wed May 18 14:18:39 2017

@author: mullerrs

RUN SCRIPT IN SAME DIRECTORY AS ALL .dat DOCUMENTS OF THE DOM MEASUREMENTS
INCLUDING THE ELECTRONICS .wav FILE OF THE NOISE MEASUREMENT, AND
INCLUDING THE REFERENCE HYDROPHONE .dat MEASUREMENT TO CALIBRATE WITH
I'm sorry I didn't made this automatic yet.

WHEN RUNNING SCRIPT IN COMMANDLINE IT DOES NOT NEED ANY ARGUMENTS (yet)

I also want to implement a -h function and other parsers... 
but first things first.

"""
#==============================================================================
# IMPORT LIBRARIES, FUNCTIONS AND CLASSES
#==============================================================================

import sys
import numpy as np
import struct
import os
import math
import scipy.io
import matplotlib.pyplot as plt
from datetime import datetime
from scipy.signal import butter, lfilter
import wavio
from numpy.linalg import norm as Norm
from numpy.fft import fft

#==============================================================================
# GLOBAL VARIABLES
#==============================================================================

audioFs     = 15.25 # frequency of audio source pulses

save_txt    = False # if save_txt == True --> all txt files, are being saved
plotjes     = True # if plotjes == True --> plots are being made                <== must be true to save png
save_png    = False # if save_png == True --> all png files, are being saved
show        = True # if show == True --> all plot's are being showed

#==============================================================================
# FILTER
#==============================================================================

def butter_lowpass(cutoff, Fs, order=5):
    nyq = 0.5 * Fs
    normal_cutoff = cutoff / nyq
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    return b, a

def butter_lowpass_filter(data, basename, cutoff, Fs, order=5):
    b, a = butter_lowpass(cutoff, Fs, order=order)
    y = lfilter(b, a, data)
    basename = basename+'_lp%s' %(cutoff/1000)
    return y, basename

def butter_highpass(cutoff, Fs, order=5):
    nyq = 0.5 * Fs
    normal_cutoff = cutoff / nyq
    b, a = butter(order, normal_cutoff, btype='high', analog=False)
    return b, a

def butter_highpass_filter(data, basename, cutoff, Fs, order=5):
    b, a = butter_highpass(cutoff, Fs, order=order)
    y = lfilter(b, a, data)
    basename = basename+'_hp%s' %(cutoff/1000)
    return y, basename

def butter_bandpass_filter(data, basename, lowcut, highcut, Fs, order=5):
    y = butter_lowpass_filter(data, basename, highcut, Fs, order)[0]
    basename = basename+'_bp%s-%s' %(lowcut/1000, highcut/1000)
    return butter_highpass_filter(y, basename, lowcut, Fs, order)[0], basename

#==============================================================================
# CUT
#==============================================================================

def Cut(data, samples, Fs, basename):                                           # <== idealiter overal zelfdes soort volgorde aanhouden

    MAX = np.max(data)                 # Maximum valule in dataset
    treshold = MAX/2.5
    WHERE = np.where(data > treshold)  # Samples with value above treshold
    WHERE = np.asarray(WHERE[0])       # Make an array of the tuple
    start_pulse = [WHERE[0]]           # Create an array start_pulse & put 
                                       # first sample above MAX/2 init
    # In de komende loop: check voor alle samples boven MAX/2 of ze bij 
    # dezelfde pulse horen, of een nieuwe pulse zijn definitie nieuwe pulse is 
    # als twee opeenvolgende samples met waarde boven treshold, meer dan Fs/16 
    # van elkaar af liggen als een nieuwe pulse gevonden is: voeg deze toe aan 
    # start_pulse
    i = 1
    while i < len(WHERE):
        if WHERE[i] > (start_pulse[-1] + (Fs/(int(audioFs)+1))):
            start_pulse.append(WHERE[i])
        i += 1

    N = samples
    cutted_data = np.array([])
    i = 0
    j = 0 # <- to count the removed pulses
    while i < len(start_pulse):
        range_I = start_pulse[i] - N/2
        cutted_pulse = (data[range_I : range_I+N])    #NB: hij heeft links 1 
        # sample meer dan rechts van start_pulse. : is TOT, niet t/m.

        if range_I < 0:
            print 'IndexError occurred - Nothing to worry about: Skipped pulse'
            i += 1
            j += 1
        elif range_I+N > len(data):
            print 'IndexError occurred - Nothing to worry about: Skipped pulse'
            i += 1
            j += 1
        else:
            cutted_data = np.concatenate((cutted_data, cutted_pulse))
            i += 1

    basename = basename + '_pulsecut%s' %samples

    x = float(len(cutted_data))/float(N) #how often 1024 samples in cutted data
    if x != len(start_pulse) - j:
        print 'Not all pulses are cutted correctly'
        print '\n Extra info to search for error: \n'
        print 'Expected amount of pulses from data: int of:      ', len(data)/Fs*audioFs
        print 'Amount of pulses that should have been plotted =  ', len(start_pulse)
        print 'Amount of pulses cutted =                         ', x
        print 'Max value                                         ', MAX
        print 'Treshold value                                    ', MAX/treshold
        print 'Start of the pulses that should have been plotted: \n', start_pulse

    ####### PLOT DATA #######
    if plotjes == True:
        l = N      # unit of file to read in:  [raw: Fs],   [cut: 1024]
        x = np.linspace(0, l-1, l)
        plt.plot(x, cutted_data[0:l], linewidth=1)
        plt.tick_params(axis='both', labelsize=15)
        if '_wavelet_' in basename:
            HYD = 'DOM piezo'
        else:
            HYD = 'reference hydrophone'
        plt.title('\n Single pulse of %s' %HYD, fontsize=24)
        plt.suptitle('used: %s'%basename, fontsize=9)
        plt.xlabel('# samples', fontsize=18)
        plt.ylabel('amp', fontsize=18)
        plt.tight_layout()
        plt.grid()
        plt.grid('on','minor')
        if save_png == True:
            plt.savefig(basename+'_data_one-pulse.png') #_data_1wav_zoom.png
        if show == True:
            plt.show()
        plt.clf()
    #########################

    return cutted_data, basename

#==============================================================================
# Make time span
#==============================================================================

def MakeTimeSpan(length, Fs):
    return np.linspace(0, length-1, length)/Fs

#==============================================================================
# Make frequency span
#==============================================================================

def MakeFreqSpan(Fs, NFFT):
    return np.linspace(0, NFFT/2-1, NFFT/2)*float(Fs)/float(NFFT)

#==============================================================================
# Plot Ed
#==============================================================================
def plot(data, basename, Fs, NFFT):
    length  = np.size(data)
    length  = len(data)
    time    = MakeTimeSpan(length, Fs)

    plt.plot(time, data)
    stdev = (np.std(data))
    str = "$\sigma$: %6.4f %s" % (stdev, 'Pa');
    plt.tick_params(axis='both', labelsize=15)
    plt.title('\n\nTimetrace', fontsize=24)
    plt.suptitle('used: %s\nNFFT: %s\n \%s'%(basename, NFFT, str), fontsize=9)
    str = "%s -> [%s] " %('amplitude', 'Pa');
    plt.ylabel(str, fontsize=18)
    plt.autoscale(tight = True)
    plt.xlabel('time ->[s]', fontsize=18)
    plt.grid()
    plt.grid('on', 'minor')
    plt.tight_layout()
    if save_png == True:
        plt.savefig(basename+'.png')
    if show == True:
        plt.show()
    plt.clf()

########################                                                        #<== is nu een beetje dubbelop

def MakeSpectrum(Fs, NFFT, data, basename):

    length     = len(data)
    PSD        = np.zeros((NFFT/2, 1), dtype=float)
    freq       = MakeFreqSpan(Fs, NFFT)
    Segment    = int(length/NFFT)
    wnd        = np.hanning(NFFT)
    norm       = Norm(wnd)**2
    double2single   = 2.0
    for span in range(0, Segment):
        Bg                 = span*NFFT
        end                = Bg+NFFT
        yw                 = wnd*data[Bg:end]
        a                  = fft(yw, NFFT)
        ac                 = np.conj(a)
        pxx                = np.abs(a*ac)
        PSD[:, 0]         +=  double2single*pxx[0:NFFT/2]
    PSD[:, 0]  /= (float(Segment)*NFFT*norm)                                    #<== this method of averaging several spectra is also known under names such as "Welch's method ... [Hei02]
    basename = basename+ '_spectrum'
    return PSD[:, 0], basename

########################

def PlotSpectrum(Fs, NFFT, data, basename):
#    if "10-60kHz" in basename:
#        FcumBegin   = 10e3 - 5e3
#        FcumEnd     = 60e3 + 5e3
#    elif "50_80kHz" in basename:
#        FcumBegin   = 50e3 - 5e3
#        FcumEnd     = 80e3 + 5e3
#    elif "Sw50to80kHz" in basename:
#        FcumBegin   = 50e3 - 5e3
#        FcumEnd     = 80e3 + 5e3
#    elif "DOM_40_ruis" in basename:
#        FcumBegin   = 30e3 - 5e3
#        FcumEnd     = 30e3 + 5e3
#    elif "2kHz" in basename:
#        FcumBegin   = 2e3 - 5e3
#        FcumEnd     = 2e3 + 5e3
#    elif "_10kHz" in basename:
##        FcumBegin   = 10e3 - 5e3
##        FcumEnd     = 10e3 + 5e3
#        FcumBegin   = 3500
#        FcumEnd     = 22500
#    elif "_20kHz" in basename:
#        FcumBegin   = 20e3 - 5e3
#        FcumEnd     = 20e3 + 5e3
#    elif "_30kHz" in basename:
#        FcumBegin   = 30e3 - 5e3
#        FcumEnd     = 30e3 + 5e3
#    elif "_40kHz" in basename:
#        FcumBegin   = 40e3 - 5e3
#        FcumEnd     = 40e3 + 5e3
#    elif "_50kHz" in basename:
#        FcumBegin   = 50e3 - 5e3
#        FcumEnd     = 50e3 + 5e3
#    elif "_60kHz" in basename:
#        FcumBegin   = 60e3 - 5e3
#        FcumEnd     = 60e3 + 5e3
#    else:
#        FcumBegin   = 3500
#        FcumEnd     = Fs/2
    FcumBegin = 5000
    FcumEnd = 80000

    PSD      = data    
    dF       = float(Fs)/float(NFFT)    # amount of freq per bin
    idxB     = int(FcumBegin/dF)        # bin to start cumsum
    idxE     = int(FcumEnd/dF)          # bin to stop cumsum
    cummulative = PSD[idxB:idxE].cumsum();
    std = 10*np.log10(cummulative[-1]); #According to Parseval
    freq = MakeFreqSpan(Fs, NFFT)

#################################################
#    Plot sqrt(PSD), dus met Pa op y-as
#################################################
#    sqrt_std = np.sqrt(cummulative[-1])
#    plt.plot(freq, np.sqrt(abs(PSD)), zorder = 1, linewidth=1) # PSD in dB re amp^2
#    #plot cummulatives
#    plt.plot(freq[idxB:idxE], abs(cummulative), linewidth=1)
#    plt.plot(freq[idxE:idxB-1:-1], abs(PSD[idxE:idxB-1:-1].cumsum()), linewidth=1)
#    plt.tick_params(axis='both', labelsize=15)
#    str = "$\sigma$ = %6.4f Pa" % (sqrt_std)
#    plt.title('\n\nPowerspectrum', fontsize=24)
#    plt.suptitle('used: %s\nNFFT: %s\n %s' %(basename, NFFT, str), fontsize=9)
#    plt.ylabel("ampl [Pa] ", fontsize=18)
#    plt.autoscale(tight = True)
#    plt.grid()
#    plt.grid('on', 'minor')
#    plt.xlabel('Frequency [Hz]', fontsize=18)
#    plt.tight_layout()
#    if save_png == True:
#        plt.savefig(basename+'.png')
#    if show == True:
#        plt.show()
#    plt.clf()

#################################################
#   Plot 10*log10(PSD), dus dB re Pa op y-as
#################################################
    plt.plot(freq, 10.*np.log10(abs(PSD)), linewidth=1) # PSD in dB re amp^2
    #plot cummulatives
    plt.plot(freq[idxB:idxE],10.*np.log10(abs(cummulative)), linewidth=1)
    plt.plot(freq[idxE:idxB-1:-1], 10.*np.log10(abs(PSD[idxE:idxB-1:-1].cumsum())), linewidth=1)
    plt.tick_params(axis='both', labelsize=15)
    str = "$\sigma$ = %6.4f dB re Pa" % (std)
    plt.title('\n\nPowerspectrum', fontsize=24)
    plt.suptitle('used: %s\nNFFT: %s\n %s' %(basename, NFFT, str), fontsize=9)
    plt.ylabel("dB re Pa", fontsize=18)
    plt.autoscale(tight = True)
    plt.grid()
    plt.grid('on', 'minor')
    plt.xlabel('Frequency [Hz]', fontsize=18)
    plt.tight_layout()
    if save_png == True:
        plt.savefig(basename+'_dBrePa.png')
    if show == True:
        plt.show()
    plt.clf()

#################################################
#   Plot 10*log10(PSD), dus dB re microPa op y-as
#################################################

    std = 10*np.log10(cummulative[-1]/(1e-6)**2);
    plt.plot(freq, 20*np.log10(np.sqrt(abs(PSD))/1e-6), linewidth=1) # dB re $\mu$Pa
    #plot cummulatives
    plt.plot(freq[idxB:idxE],10.*np.log10(abs(cummulative)/(1e-6)**2), linewidth=1)
    plt.plot(freq[idxE:idxB-1:-1], 10.*np.log10(abs(PSD[idxE:idxB-1:-1].cumsum())/(1e-6)**2), linewidth=1)
    plt.tick_params(axis='both', labelsize=15)
    str = "$\sigma$ = %6.4f dB re Pa" % (std)
    plt.title('\n\nPowerspectrum', fontsize=24)
    plt.suptitle('used: %s\nNFFT: %s\n %s' %(basename, NFFT, str), fontsize=9)
    plt.ylabel("dB re $\mu$Pa", fontsize=18)
    plt.autoscale(tight = True)
    plt.grid()
    plt.grid('on', 'minor')
    plt.xlabel('Frequency [Hz]', fontsize=18)
    plt.tight_layout()
    if save_png == True:
        plt.savefig(basename+'_dBremuPa.png')
    if show == True:
        plt.show()
    plt.clf()

    return

#==============================================================================
# PSD Ed - for calibration
#==============================================================================

def Spectrum(Fs, NFFT, data, basename):

    length     = len(data)
    PSD        = np.zeros((NFFT/2, 1), dtype=float)
    freq       = MakeFreqSpan(Fs, NFFT)
    Segment    = int(length/NFFT)
    wnd        = np.ones(NFFT) 
    wnd        = np.hanning(NFFT)
    norm       = Norm(wnd)**2
    double2single   = 2.0
    for span in range(0, Segment):
        Bg                 = span*NFFT
        end                = Bg+NFFT
        yw                 = wnd*data[Bg:end]
        a                  = fft(yw, NFFT)
        ac                 = np.conj(a)
        pxx                = np.abs(a*ac)
        PSD[:, 0]         +=  double2single*pxx[0:NFFT/2]
    PSD[:, 0]  /= (float(Segment)*NFFT*norm)

    if plotjes == True:

        plt.plot(freq, PSD[:,0], linewidth=1)
        plt.tick_params(axis='both', labelsize=15)
        plt.title('\nPowerspectrum', fontsize=24)
        plt.suptitle('used: %s\nNFFT: %s'%(basename, NFFT), fontsize=9)
        plt.xlabel('Frequency [Hz]', fontsize=18)
        plt.ylabel('amp$^{2}$', fontsize=18)
        plt.tight_layout()
        plt.grid()
        plt.grid('on','minor')
        if save_png == True:
            plt.savefig(basename+'_amp2_calibration.png')
        if show == True:
            plt.show()
        plt.clf()

        plt.plot(freq, np.sqrt(PSD[:,0]), linewidth=1)
        plt.tick_params(axis='both', labelsize=15)
        plt.title('\nPowerspectrum', fontsize=24)
        plt.suptitle('used: %s\nNFFT: %s'%(basename, NFFT), fontsize=9)
        plt.xlabel('Frequency [Hz]', fontsize=18)
        plt.ylabel('amp', fontsize=18)
        plt.tight_layout()
        plt.grid()
        plt.grid('on','minor')
        if save_png == True:
            plt.savefig(basename+'_amp_calibration.png')
        if show == True:
            plt.show()
        plt.clf()

        plt.plot(freq, 20*np.log10(np.sqrt(PSD[:,0])), linewidth=1)
        plt.tick_params(axis='both', labelsize=15)
        plt.title('\nPowerspectrum', fontsize=24)
        plt.suptitle('used: %s\nNFFT: %s'%(basename, NFFT), fontsize=9)
        plt.xlabel('Frequency [Hz]', fontsize=18)
        plt.ylabel('dB re amp', fontsize=18)
        plt.tight_layout()
        plt.grid()
        plt.grid('on','minor')
        if save_png == True:
            plt.savefig(basename+'_dBreamp_calibration.png')
        if show == True:
            plt.show()
        plt.clf()

        plt.plot(freq, 20*np.log10(np.sqrt(PSD[:,0])/(1e-6)), linewidth=1)
        plt.tick_params(axis='both', labelsize=15)
        plt.title('\nPowerspectrum', fontsize=24)
        plt.suptitle('used: %s\nNFFT: %s'%(basename, NFFT), fontsize=9)
        plt.xlabel('Frequency [Hz]', fontsize=18)
        plt.ylabel('dB re $\mu$amp', fontsize=18)
        plt.tight_layout()
        plt.grid()
        plt.grid('on','minor')
        if save_png == True:
            plt.savefig(basename+'_dBremuamp_calibration.png')
        if show == True:
            plt.show()
        plt.clf()

    basename = basename + '_spectrum'                                           # <== nog zeggen wat voor een type spectrum? db / amp2 / amp

    return freq, PSD[:,0], basename

#==============================================================================
# GET RATIO ELECTRONICS                                                         <==
#==============================================================================

def Elec_ratio(wav_file, cal_freq, NFFT, lowcut, highcut):

    basename = os.path.splitext(os.path.basename(wav_file))[0]
    print('- Electronics of the piezo is being examined using:\n   %s' %basename)

    # get properties of the file
    data                = wav_file
    wavo                = wavio.read(data)  # returns data.shape, dtype, rate 
                                            # and sampwitdt (sampwidth == 3
                                            # <-> 24 bits WAV)
    data                = wavo.data[:,0]    # this can be saved as a .dat file
    Fs                  = wavo.rate
    length              = wavo.data.shape[0]
    NFFT                = NFFT

    # filter
    data = butter_bandpass_filter(data, basename, lowcut, highcut, Fs, order=5)[0]

    # Make - PSD - amp2
    freq, PSD, basename = Spectrum(Fs, NFFT, data, basename)

    # Normalize
    sqrt_PSD = np.sqrt(PSD)
    PSD10kHz = sqrt_PSD[NFFT*cal_freq/Fs]    # wat is de waarde bij cal_freq?
    ratios = sqrt_PSD/PSD10kHz
    basename = basename + '_normalized'

    # save
    if save_txt == True:    
        np.savetxt('20170419_6_white_noise_0.9Vrms_RATIO_nomalized.dat', ratios)

    # plot ratios
    if plotjes == True:
        plt.plot(freq, ratios, linewidth=1)
        plt.tick_params(axis='both', labelsize=15)
        plt.title('\n\n\nNormalized ratio electronics', fontsize=24)
        plt.suptitle('Freq at calibration freq %skHz equals one' %int(cal_freq/1000.))
        plt.xlabel('Frequency [Hz]', fontsize=18)
        plt.ylabel('amp', fontsize=18) # only /$\sqrt(Hz)$ if PDS so NFFT = Fs
        plt.ylim(0, 10)
        plt.grid()
        plt.grid('on','minor')
        plt.tight_layout()
        if save_png == True:
            plt.savefig(basename+'_amp_calibration.png')
        if show == True:
            plt.show()
        plt.clf()

    # plot ratios on dB scale
    if plotjes == True:
        plt.plot(freq, 20*np.log10(ratios), linewidth=1)
        plt.tick_params(axis='both', labelsize=15)
        plt.title('\n\n\nNormalized ratio electronics')
        plt.suptitle('Freq at calibration freq %skHz equals one' %int(cal_freq/1000.))
        plt.xlabel('Frequency [Hz]', fontsize=18)
        plt.ylabel('dB re amp',fontsize=18) # only sqrt(Hz) if PDS so NFFT = Fs
        plt.ylim(-20, 20)
        plt.grid()
        plt.grid('on','minor')
        plt.tight_layout()
        if save_png == True:
            plt.savefig(basename+'_dBreamp_calibration.png')
        if show == True:
            plt.show()
        plt.clf()

    return freq, ratios

#==============================================================================
# CORRECTION FOR DISTANCE
#==============================================================================

def distance_correction(config):

    print('- Determine distance correction\n'
          '   based on 1 over distance squared relation')
    
    if config == 1:
        RD = 0          #[m]
        SD = 1          #[m]
        SR = SD - RD    #[m]
        print 'distance_correction has not been determined for this config'
    elif config == 2:
        RD = 0.52       #[m]
        SD = 0.97       #[m]
        SR = SD - RD    #[m]
    elif config == 3:
        RD = 0          #[m]
        SD = 1          #[m]
        SR = SD - RD    #[m]
        print 'distance_correction has not been determined for this config'
    elif config == 4:
        RD = 0          #[m]
        SD = 1          #[m]
        SR = SD - RD    #[m]
        print 'distance_correction has not been determined for this config'
    elif config == 5:
        RD = 0          #[m]
        SD = 1          #[m]
        SR = SD - RD    #[m]
        print 'distance_correction has not been determined for this config'
    elif config == 6:
        RD = 0          #[m]
        SD = 1          #[m]
        SR = SD - RD    #[m]
        print 'distance_correction has not been determined for this config'

    distance_correction = SR**2/SD**2

    return distance_correction

#==============================================================================
# CALIBRATION VALUE                                                             <==
#==============================================================================

def IJk_value(hyd_data, basename, Fs, NFFT, FcumBegin, FcumEnd):

    Fs = Fs

    if Fs == int(25e6/128):
        print '- Determine calibration value of the DOM using:\n   %s' %basename
    elif Fs == 2e5:
        print '- Determine calibration value of the REF using:\n   %s' %basename

    # PSD in amp2
    freq, PSD, basename = Spectrum(Fs, NFFT, hyd_data, basename)

    # Because noise measurement only shows peaks of electronics
    # we do not subtract noise.

    # determine cumsum
    dF       = float(Fs)/float(NFFT)    # amount of freq per bin
    idxB     = int(FcumBegin/dF)        # bin to start cumsum
    idxE     = int(FcumEnd/dF)          # bin to stop cumsum
    cummulative = PSD[idxB:idxE].cumsum();
    std = 10*np.log10(cummulative[-1]); #According to Parseval
    sqrt_std = np.sqrt(cummulative[-1])

    print '   sigma = %6.4f dB re amp, or sigma = %6.4f amp' %(std, sqrt_std)
    
#WAT ALS IK GEEN ABS DOE?
    if plotjes == True:

        plt.plot(freq, np.sqrt(abs(PSD)), zorder = 1, linewidth=1) # PSD in dB re amp^2
        #plot cummulatives
        plt.plot(freq[idxB:idxE], abs(cummulative), linewidth=1)
        plt.plot(freq[idxE:idxB-1:-1], abs(PSD[idxE:idxB-1:-1].cumsum()), linewidth=1)
        plt.tick_params(axis='both', labelsize=15)
        str = "$\sigma$ = %6.4f amp" % (sqrt_std)
        plt.title('\n\nPowerspectrum', fontsize=24)
        plt.suptitle('used: %s\nNFFT: %s\n %s' %(basename, NFFT, str), fontsize=9)
        plt.ylabel("amp ", fontsize=18)
        plt.autoscale(tight = True)
        plt.grid()
        plt.grid('on', 'minor')
        plt.xlabel('Frequency [Hz]', fontsize=18)
        plt.tight_layout()
        if save_png == True:
            plt.savefig(basename+'.png')
        if show == True:
            plt.show()
        plt.clf()

        plt.plot(freq, 10.*np.log10(abs(PSD)), zorder = 1, linewidth=1) # PSD in dB re amp^2
        #plot cummulatives
        plt.plot(freq[idxB:idxE],10.*np.log10(abs(cummulative)), linewidth=1)
        plt.plot(freq[idxE:idxB-1:-1], 10.*np.log10(abs(PSD[idxE:idxB-1:-1].cumsum())), linewidth=1)
        plt.tick_params(axis='both', labelsize=15)
        str = "$\sigma$ = %6.4f dB re amp" % (std)
        plt.title('\n\nPowerspectrum', fontsize=24)
        plt.suptitle('used: %s\nNFFT: %s\n %s' %(basename, NFFT, str), fontsize=9)
        plt.ylabel("dB re amp", fontsize=18)
        plt.autoscale(tight = True)
        plt.grid()
        plt.grid('on', 'minor')
        plt.xlabel('Frequency [Hz]', fontsize=18)
        plt.tight_layout()
        if save_png == True:
            plt.savefig(basename+'_dB.png')
        if show == True:
            plt.show()
        plt.clf()

    calibration_value = sqrt_std
    
    return calibration_value

#==============================================================================
# PREPARE THE FILES USED FOR CALIBRATIONVALUES
#==============================================================================

def get_cut_filt_IJkfile(file_ijkHYD, Fs, lowcut, highcut, NFFT):

    Fs = Fs
    ijkHYD = np.genfromtxt(file_ijkHYD)
    basename = os.path.splitext(os.path.basename(file_ijkHYD))[0]
    ijkHYD, basename = butter_bandpass_filter(ijkHYD, basename, lowcut, highcut, Fs, order=5)
    ijkHYD, basename = Cut(ijkHYD, NFFT, Fs, basename)

    return ijkHYD, basename
    
#==============================================================================
# CALIBRATE
#==============================================================================

def CalibrationArray(file_electronics, file_ijkDOM, file_ijkREF, lowcut, highcut, cal_freq, NFFT, Fs_DOM, Fs_REF, cal_min, cal_max):
#                       .wav file       .dat file    .dat file     [Hz]    [Hz]     [Hz]  [samples]
    print('\nGOING TO MAKE CALIBRATION ARRAY')

    elec_freq, elec_ratio = Elec_ratio(file_electronics, cal_freq, NFFT, lowcut, highcut)
#   dimension less, since it's a ratio

    ijkDOM, basenameDOM = get_cut_filt_IJkfile(file_ijkDOM, Fs_DOM, lowcut, highcut, NFFT)
    ijkREF, basenameREF = get_cut_filt_IJkfile(file_ijkREF, Fs_REF, lowcut, highcut, NFFT)

    DOM_val_X = IJk_value(ijkDOM, basenameDOM, int(25e6/128), NFFT, cal_min, cal_max) # in [X]
    REF_val_Pa = IJk_value(ijkREF, basenameREF, int(2e5), NFFT, cal_min, cal_max) # in [Pa]

    sensitivity_value = (DOM_val_X / REF_val_Pa)   # in [X]/[Pa]
    dist_cor = distance_correction(config = 2)
    calibration_value = sensitivity_value / dist_cor

    print '- Used values used for calibration:'
    print '   Sensitivity value (DOM/REF) = %s' %sensitivity_value
    print '   Distance correction value   = %s' %dist_cor
    print '   Resulting calibration value = %s' %calibration_value

    return elec_freq, elec_ratio, calibration_value                                        # Wil je niet liever x en y?

#==============================================================================
# MAIN
#==============================================================================

def main():
    #start timer
    startTime = datetime.now()

    # input values for calubration:
    cal_type        = 'wav10kHz'                                                # Specify type of calibration for savingname.
    file_electronics= '20170419_6_white_noise_0.9Vrms.WAV'           # .wav     # 
    file_ijkDOM     = '20170310_config2_DOM_6_v1_wavelet_10kHz.dat'  # .dat     # 
    file_ijkREF     = 'Time_Conf2Domwavelet10kHzv1.dat'              # .dat     # 
    lowcut          = 1000                                           # [Hz]     # 
    highcut         = 80000                                          # [Hz]     # 
    cal_freq        = 10000                                          # [Hz]     # 
    samples_to_cut  = 1024                                           # [samp]   # 
    NFFT            = samples_to_cut # do NOT change this!           # [samp]   # 
    Fs_DOM          = int(25e6/128)                                  # [samp/s] # 
    Fs_REF          = int(2e5)                                       # [samp/s] # 
    cal_min         = 3500                                           # [Hz]     # 
    cal_max         = 22500                                          # [Hz]     # 

    #get electronics ratio, and calibration value
    #elec_freq, elec_ratio, calibration_value = CalibrationArray(file_electronics, file_ijkDOM, file_ijkREF, lowcut, highcut, cal_freq, NFFT, Fs_DOM, Fs_REF, cal_min, cal_max)
    #due to different Fs, elec_freq is slightly different from freq

    # DOM files
    files   = ['20170309_config1_DOM_7_toon_30kHz.dat', '20170309_config1_DOM_10_burst_30kHz.dat',  '20170310_config2_DOM_10_v1_wavelet_30kHz.dat']

    # all Rectangular window!    
#    files   = ['sine', 'Burst_30kHz_5Vpp_12cycles_hamming_AWG33220.dat', 'Swp_5-65kHz_5Vpp_Clock2MHz_windowed_AWG33220.dat', '20170309_config2_DOM_17_sweep10-60kHz_10Vpp.dat', 'Wavelet_30kHz_5Vpp.dat']
    typ     = ['Tone 10 kHz', 'Burst 10kHz', 'Sweep 10-60kHz', 'Wavelet 10kHz']

    for i in files:
        name = i
        if i =='sine':
            Fs = 2e6
            F0 = 30000
            length = 1024
            time = np.linspace(0, length-1, length)/Fs
            data = 6e5*np.sin(2*np.pi*F0*time)
        elif i == '20170309_config1_DOM_7_toon_30kHz.dat':
            data = np.genfromtxt(i, max_rows = 2e4) # <== don't read in the whole file --> gives OverflowError
        else:
            data = np.genfromtxt(i)

        length = len(data)

        if 'sine' in i:
            filename = 'Tone' #'Tone 30kHz'
            Fs = Fs
        elif 'toon' in i:
            filename = 'Tone' #'Tone 30kHz'
            Fs = int(25e6/128)
        elif 'Burst' in i:
            filename = 'Burst'# 'Burst 30kHz'
            Fs = 2e6
        elif 'burst' in i:
            filename = 'Burst'# 'Burst 30kHz'
            Fs = int(25e6/128)
        elif 'Swp' in i:
            filename = 'Sweep'#'Sweep 5-65kHz'
            Fs = 2e6
        elif 'sweep' in i:
            filename = 'Sweep' #'Sweep 10-60 kHz'
            Fs = int(25e6/128)
        elif 'wavelet' in i:
            filename = 'Wavelet'#'Wavelet 30kHz'
            Fs = int(25e6/128)
        elif 'Wavelet' in i:
            filename = 'Wavelet' #'Wavelet 30kHz'
            Fs = 1e6
        else:
            print 'not sure how to handle file...'

#==============================================================================
# SIMPLE
#==============================================================================
        time = np.linspace(0, length-1, length)/Fs
        if 'Swp' in i:        
            plt.plot(time[0:1200]/(1e-3), data[0:1200]/(1e5), linewidth=1)
        else:
            plt.plot(time/(1e-3), data/(1e5), linewidth=1)
        plt.title('\n%s' %filename, fontsize=24)
        plt.tick_params(axis='both', labelsize=15)
        plt.xlabel('Time [ms]', fontsize=18)
        plt.yticks
        plt.autoscale(tight = True)
        plt.grid()
        plt.grid('on', 'minor')
        plt.tight_layout()
        plt.show()
        plt.clf()

        freq = MakeFreqSpan(Fs, NFFT)
        if 'Swp' in i:        
            PSD = Spectrum(Fs, NFFT, data[0:1200], 'basename')[1]
        else:
            PSD = Spectrum(Fs, NFFT, data, 'basename')[1]
        plt.plot(freq, 20*np.log10(np.sqrt(PSD)), linewidth=1)
        plt.title('\n%s' %filename, fontsize=24)
        plt.tick_params(axis='both', labelsize=15)
        plt.xlabel('Frequency [Hz]', fontsize=18)
        plt.autoscale(tight = True)
        plt.grid()
        plt.grid('on', 'minor')
        plt.tight_layout()
        plt.show()
        plt.clf()

#==============================================================================
#  MORE TEXT
#==============================================================================

#        time = np.linspace(0, length-1, length)/Fs
#        plt.plot(time/(1e-3), data/(1e5), linewidth=1)
#        plt.title('\n\nTimetrace \n%s' %filename, fontsize=24)
#        plt.suptitle('used: %s' %(name), fontsize=9)
#        plt.tick_params(axis='both', labelsize=15)
#        plt.xlabel('Time [ms]', fontsize=18)
#        plt.ylabel('Amplitude [x 10$^{5}$]', fontsize=18)
#        plt.autoscale(tight = True)
#        plt.grid()
#        plt.grid('on', 'minor')
#        plt.tight_layout()
#        plt.show()
#        plt.clf()
#
#        freq = MakeFreqSpan(Fs, NFFT)
#        PSD = Spectrum(Fs, NFFT, data, 'basename')[1]
#        
#        plt.plot(freq, 20*np.log10(np.sqrt(PSD)), linewidth=1)
#        plt.title('\n\nSpectrum \n%s' %filename, fontsize=24)
#        plt.suptitle('used: %s' %(name), fontsize=9)
#        plt.tick_params(axis='both', labelsize=15)
#        plt.xlabel('Time [s]', fontsize=18)
#        plt.ylabel('dB re amplitude', fontsize=18)
#        plt.autoscale(tight = True)
#        plt.grid()
#        plt.grid('on', 'minor')
#        plt.tight_layout()
#        plt.show()
#        plt.clf()

    raise sys.exit()

    #end timer
    print datetime.now() - startTime

if __name__ == "__main__":
    sys.exit(main())
