# -*- coding: utf-8 -*-
"""
                                                                                Created on Wed May 18 14:18:39 2017

@author: mullerrs

RUN SCRIPT IN SAME DIRECTORY AS ALL .dat DOCUMENTS OF THE DOM MEASUREMENTS
INCLUDING THE ELECTRONICS .wav FILE OF THE NOISE MEASUREMENT, AND
INCLUDING THE REFERENCE HYDROPHONE .dat MEASUREMENT TO CALIBRATE WITH

WHEN RUNNING SCRIPT IN COMMANDLINE IT DOES NOT NEED ANY ARGUMENTS
                                                                                I do want to implement a -h function and other parsers... 
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

audioFs     = 15.25 # frequency of audio source pulses [Hz]

save_txt    = False # if save_txt == True   --> all txt files, are being saved
plots       = True # if plots == True      --> plots are being made
# plots must be true to be able to use next 3 options:
save_png    = False # if save_png == True   --> all png files, are being saved
show        = True # if show == True       --> all plot's are being showed
cumm        = False # if cumm == Trus       --> cumpute cummulatives.               <== does not always work (like in Ed's script)

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
    print basename
    return butter_highpass_filter(y, basename, lowcut, Fs, order)[0], basename

#==============================================================================
# CUT
#==============================================================================

def Cut(data, samples, Fs, basename):                                           # <== idealiter overal zelfdes soort volgorde aanhouden

    MAX = np.max(data)                 # Maximum valule in dataset
    treshold = MAX/2.5
    WHERE = np.where(data > treshold)  # Samples with value above treshold
    WHERE = np.asarray(WHERE[0])       # Make an array of the tuple
    start_pulse = [WHERE[0]]
    # start_pulse now is an array with every first element of a new pulse.
    # 'first element' as first time above defined treshold.

    # In the next loop: check for all samples above the treshold if they belong
    # to the same pulse, or if they belong to a new pulse. Definition of new is
    # when two following samples above the treshold  are separated by more than 
    # Fs/16 samples. When a new pulse is being found: add it to the array:
    # start_pulse
    i = 1
    while i < len(WHERE):
        if WHERE[i] > (start_pulse[-1] + (Fs/(int(audioFs)+1))):
            start_pulse.append(WHERE[i])
        i += 1

    N = Fs/50
    cutted_data = np.array([])
    i = 0                               # to loop through start_pulse
    j = 0                               # to count the removed pulses
    while i < len(start_pulse):
        range_I = start_pulse[i] - N/2  # Strat of pulse should be in middle
        cutted_pulse = (data[range_I : range_I+N])    
        # This will result in 1 sample more on the right side, compared to left

    # In the next if-loop, pulses will be removed, if the N samples around the
    # first or last pulse fall outside the measured data. Otherwhise the 
    # cutted_pulse is added to the cutted data.
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

    # This was to check the cutting procedure
    x = float(len(cutted_data))/float(N) # Amount cutted pulses there should be
    if x != len(start_pulse) - j:        # Check if this is indeed true
        print 'Not all pulses are cutted correctly'
        print '\n Extra info to search for error: \n'
        print ('Expected amount of pulses from data: int of:      ',
               len(data)/Fs*audioFs)
        print ('Amount of pulses that should have been plotted =  ',
               len(start_pulse))
        print ('Amount of pulses cutted =                         ',
               x)
        print ('Max value                                         ',
               MAX)
        print ('Treshold value                                    ',
               MAX/treshold)
        print ('Start of the pulses that should have been plotted: \n',
               start_pulse)

#    ####### PLOT DATA #######
#    # amount of samples to plot is defined with l:
#    if plots == True:
#        l = N      # amount of samples to read in                              <== [raw: Fs],   [cut: 1024]
#        x = np.linspace(0, l-1, l)
#        plt.plot(x, cutted_data[0:l], linewidth=1)
#        plt.tick_params(axis='both', labelsize=15)
#        # define the hydrophone that was used: DOM or Reference hydrophone        
#        # assumes only wavelet is being cutted.
#        if '_wavelet_' in basename:
#            HYD = 'DOM piezo'
#        else:
#            HYD = 'reference hydrophone'
#        plt.title('\n Single pulse of %s' %HYD, fontsize=24)
#        plt.suptitle('used: %s'%basename, fontsize=9)
#        plt.xlabel('# samples', fontsize=18)
#        plt.ylabel('amp', fontsize=18)
#        plt.tight_layout()
#        plt.grid()
#        plt.grid('on','minor')
#        if save_png == True:
#            plt.savefig(basename+'_data_one-pulse.png')
#        if show == True:
#            plt.show()
#        plt.clf()
#    #########################

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
# based on Ed's script
#==============================================================================
def plot(data, basename, Fs, NFFT):
    length  = len(data)
    time    = MakeTimeSpan(length, Fs)

    plt.plot(time, data)
    stdev = (np.std(data))
    str = "$\sigma$: %6.4f Pa" % (stdev);
    plt.tick_params(axis='both', labelsize=15)
    plt.title('\n\nTimetrace', fontsize=24)
    plt.suptitle('used: %s\nNFFT: %s\n' %(basename, NFFT) + str, fontsize=9)
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
    if 'sweep' in basename:
        wnd    = np.ones(NFFT)
        basename = basename+'_wnd_rect'
        print '   rectangular window used'
    else:
        wnd    = np.hanning(NFFT)
        basename = basename+'_wnd_hann'
        print '   hanning window used'
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
    
    PSD      = data    
    
    if cumm == True:
        if "10-60kHz" in basename:
            FcumBegin   = 10e3 - 5e3
            FcumEnd     = 60e3 + 5e3
        elif "50_80kHz" in basename:
            FcumBegin   = 50e3 - 5e3
            FcumEnd     = 80e3 + 5e3
        elif "Sw50to80kHz" in basename:
            FcumBegin   = 50e3 - 5e3
            FcumEnd     = 80e3 + 5e3
        elif "DOM_40_ruis" in basename:
            FcumBegin   = 30e3 - 5e3
            FcumEnd     = 30e3 + 5e3
        elif "2kHz" in basename:
            FcumBegin   = 2e3 - 5e3
            FcumEnd     = 2e3 + 5e3
        elif "_10kHz" in basename:
    #        FcumBegin   = 10e3 - 5e3
    #        FcumEnd     = 10e3 + 5e3
            FcumBegin   = 3500         # better boundaries chosen by hand
            FcumEnd     = 22500        # better boundaries chosen by hand
        elif "_20kHz" in basename:
            FcumBegin   = 20e3 - 5e3
            FcumEnd     = 20e3 + 5e3
        elif "_30kHz" in basename:
            FcumBegin   = 30e3 - 5e3
            FcumEnd     = 30e3 + 5e3
        elif "_40kHz" in basename:
            FcumBegin   = 40e3 - 5e3
            FcumEnd     = 40e3 + 5e3
        elif "_50kHz" in basename:
            FcumBegin   = 50e3 - 5e3
            FcumEnd     = 50e3 + 5e3
        elif "_60kHz" in basename:
            FcumBegin   = 60e3 - 5e3
            FcumEnd     = 60e3 + 5e3
        else:
            FcumBegin   = 3500     # for lower frequencies 1/f noise dominates
            FcumEnd     = Fs/2     # Nyquist frequency

        dF       = float(Fs)/float(NFFT)    # amount of freq per bin
        idxB     = int(FcumBegin/dF)        # bin to start cummulative
        idxE     = int(FcumEnd/dF)          # bin to stop cummulative
        cummulative = PSD[idxB:idxE].cumsum(); #According to Parseval
        std = 10*np.log10(cummulative[-1]); # 10 log 10 since power spectrum
    else:
        std = ''
    freq = MakeFreqSpan(Fs, NFFT)

#################################################
#    Plot sqrt(PSD), so Pa on y-axis            # PSD in dB re amp^2
#################################################
#    sqrt_std = np.sqrt(cummulative[-1])
#    plt.plot(freq, np.sqrt(abs(PSD)), zorder = 1, linewidth=1) 
#    #plot cummulatives
#    plt.plot(freq[idxB:idxE], abs(cummulative), linewidth=1)
#    plt.plot(freq[idxE:idxB-1:-1], abs(PSD[idxE:idxB-1:-1].cumsum()), 
#        linewidth=1)
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
#   Plot 10*log10(PSD), so dB re Pa at y-axis             # PSD in dB re amp^2
#################################################
    plt.plot(freq, 10.*np.log10(abs(PSD)), linewidth=1)
    #plot cummulatives
    if cumm == True:
        plt.plot(freq[idxB:idxE], 10.*np.log10(abs(cummulative)), linewidth=1)
        plt.plot(freq[idxE:idxB-1:-1], 
                 10.*np.log10(abs(PSD[idxE:idxB-1:-1].cumsum())), linewidth=1)
    plt.tick_params(axis='both', labelsize=15)
    if cumm == True:
        str = "$\sigma$ = %6.4f dB re Pa" % (std)
    else:
        str = ''
    plt.title('\n\nPowerspectrum', fontsize=24)
    plt.suptitle('used: %s\nNFFT: %s\n %s' %(basename, NFFT, str), fontsize=9)
    plt.ylabel("dB re Pa", fontsize=18)
    #plt.autoscale('y', tight = True)
    plt.grid()
    plt.grid('on', 'minor')
    plt.xlabel('Frequency [Hz]', fontsize=18)
    plt.xlim(3000, 88000)
    plt.tight_layout()
    if save_png == True:
        plt.savefig(basename+'_dBrePa.png')
    if show == True:
        plt.show()
    plt.clf()

#################################################
#   Plot 10*log10(PSD), so dB re microPa at y-axis             # dB re $\mu$Pa
#################################################

    if cumm == True:
        std = 10*np.log10(cummulative[-1]/(1e-6)**2);
    plt.plot(freq, 20*np.log10(np.sqrt(abs(PSD))/1e-6), linewidth=1)
    #plot cummulatives
    if cumm == True:
        plt.plot(freq[idxB:idxE],10.*np.log10(abs(cummulative)/(1e-6)**2), 
                 linewidth=1)
        plt.plot(freq[idxE:idxB-1:-1], 
                 10.*np.log10(abs(PSD[idxE:idxB-1:-1].cumsum())/(1e-6)**2), 
                 linewidth=1)
    plt.tick_params(axis='both', labelsize=15)
    if cumm == True:
        str = "$\sigma$ = %6.4f dB re Pa" % (std)
    else:
        str = ''
    plt.title('\n\nPowerspectrum', fontsize=24)
    plt.suptitle('used: %s\nNFFT: %s\n %s' %(basename, NFFT, str), fontsize=9)
    plt.ylabel("dB re $\mu$Pa", fontsize=18)
    #plt.autoscale('y', tight = True)
    plt.grid()
    plt.grid('on', 'minor')
    plt.xlabel('Frequency [Hz]', fontsize=18)
    plt.xlim(3000, 88000)
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

    if plots == True:

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

    basename = basename + '_spectrum'

    return freq, PSD[:,0], basename

#==============================================================================
# GET RATIO ELECTRONICS                                                         <==
#==============================================================================

def Elec_ratio(wav_file, cal_freq, NFFT, lowcut, highcut):

    basename = os.path.splitext(os.path.basename(wav_file))[0]
    print('- Electronics of the piezo is being examined using:\n   %s' 
          %basename)

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
    data = butter_bandpass_filter(data, basename, lowcut, highcut, Fs, 
                                  order=5)[0]

    # Make - PSD - amp^2
    freq, PSD, basename = Spectrum(Fs, NFFT, data, basename)

    # Normalize
    sqrt_PSD = np.sqrt(PSD)
    PSD10kHz = sqrt_PSD[NFFT*cal_freq/Fs]    # check: value at cal_freq
    ratios = sqrt_PSD/PSD10kHz
    basename = basename + '_normalized'

    # save
    if save_txt == True:    
        np.savetxt('20170419_6_white_noise_0.9Vrms_RATIO_nomalized.dat', 
                   ratios)

    # plot ratios
    if plots == True:
        plt.plot(freq, ratios, linewidth=1)
        plt.tick_params(axis='both', labelsize=15)
        plt.title('\n\n\nNormalized ratio electronics', fontsize=24)
        plt.suptitle('Freq at calibration freq %skHz equals one' 
                      %int(cal_freq/1000.))
        plt.xlabel('Frequency [Hz]', fontsize=18)
        plt.ylabel('amp', fontsize=18) # only sqrt(Hz) if PDS so NFFT = Fs
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
    if plots == True:
        plt.plot(freq, 20*np.log10(ratios), linewidth=1)
        plt.tick_params(axis='both', labelsize=15)
        plt.title('\n\n\nNormalized ratio electronics')
        plt.suptitle('Freq at calibration freq %skHz equals one'
                      %int(cal_freq/1000.))
        plt.xlabel('Frequency [Hz]', fontsize=18)
        plt.ylabel('dB', fontsize=18) # only sqrt(Hz) if PDS so NFFT = Fs
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
    
    r_DOM = 0.44        #[m]     <==>   17 inch
    
    if config == 1:
        RD = 0.52       #[m]
        SD = 0.97       #[m]
        SR = SD - RD    #[m]
    elif config == 2:
        RD = 0.52       #[m]
        SD = 0.97       #[m]
        SR = SD - RD    #[m]
    elif config == 3:
        RD = 0.52       #[m]
        SD = 0.97       #[m]
        SR = SD - RD    #[m]
    elif config == 4:
        RD = 0.52 - r_DOM + (0.5*np.pi*r_DOM)       #[m]
        SD = 0.94 - r_DOM + (0.5*np.pi*r_DOM)       #[m]
        SR = SD - RD    #[m]
        # NB: assuming 1st pulse arrives via Glass sphere: 0.25 * 2*pi*r, 
        # not in straight path through the water
        # therefore correcting for difference in speed glass v.s. water    
    elif config == 5:
        RD = 0.52       #[m]
        SD = 0.94       #[m]
        SR = SD - RD    #[m]
        # Not sure if config is right...                                        #hing REF ongeveer tegen de DOM aan? + hoek klopt is geen 135*
        print 'distance_correction has not been determined for this config'
    elif config == 6:
        RD = 0.5 + (np.pi * r_DOM)*1500/5000. #[m]
        SD = 1.  + (np.pi * r_DOM)*1500/5000. #[m]
        SR = SD - RD    #[m]
        # NB: assuming 1st pulse arrives via Glass sphere: 0.5 * 2*pi*r, 
        # not in straight path through the inside the DOM
        # therefore correcting for difference in speed glass v.s. water
    elif congif ==7:
        RD = 1.5 + (np.pi * r_DOM)*1500/5000. #[m]
        SD = 2.  + (np.pi * r_DOM)*1500/5000. #[m]
        SR = SD - RD    #[m]
        # NB: assuming 1st pulse arrives via Glass sphere: 0.5 * 2*pi*r, 
        # not in straight path through the inside the DOM
        # therefore correcting for difference in speed glass v.s. water
    distance_correction = SR**2/SD**2

    return distance_correction

#==============================================================================
# CALIBRATION VALUE
#==============================================================================

def IJk_value(hyd_data, basename, Fs, NFFT, FcumBegin, FcumEnd):

    Fs = Fs

    if Fs == int(25e6/128):
        print ('- Determine calibration value of the DOM using:\n   %s'
                %basename)
    elif Fs == 2e5:
        print ('- Determine calibration value of the REF using:\n   %s'
                %basename)

    # PSD in amp^2
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
    
    if plots == True:
        
        # PSD in dB re amp^2
        plt.plot(freq, np.sqrt(abs(PSD)), zorder = 1, linewidth=1)
        #plot cummulatives
        plt.plot(freq[idxB:idxE], abs(cummulative), linewidth=1)
        plt.plot(freq[idxE:idxB-1:-1], abs(PSD[idxE:idxB-1:-1].cumsum()), 
                 linewidth=1)
        plt.tick_params(axis='both', labelsize=15)
        str = "$\sigma$ = %6.4f amp" % (sqrt_std)
        plt.title('\n\nPowerspectrum', fontsize=24)
        plt.suptitle('used: %s\nNFFT: %s\n %s' %(basename, NFFT, str), 
                     fontsize=9)
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

        # PSD in dB re amp^2
        plt.plot(freq, 10.*np.log10(abs(PSD)), zorder = 1, linewidth=1)
        #plot cummulatives
        plt.plot(freq[idxB:idxE],10.*np.log10(abs(cummulative)), linewidth=1)
        plt.plot(freq[idxE:idxB-1:-1], 
                 10.*np.log10(abs(PSD[idxE:idxB-1:-1].cumsum())), linewidth=1)
        plt.tick_params(axis='both', labelsize=15)
        str = "$\sigma$ = %6.4f dB re amp" % (std)
        plt.title('\n\nPowerspectrum', fontsize=24)
        plt.suptitle('used: %s\nNFFT: %s\n %s' %(basename, NFFT, str), 
                     fontsize=9)
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

    print '- Cut and filter data of:\n   %s' %file_ijkHYD

    Fs = Fs
    ijkHYD = np.genfromtxt(file_ijkHYD)
    basename = os.path.splitext(os.path.basename(file_ijkHYD))[0]
    ijkHYD, basename = butter_bandpass_filter(ijkHYD, basename, lowcut, 
                                              highcut, Fs, order=5)
    ijkHYD, basename = Cut(ijkHYD, NFFT, Fs, basename)

    return ijkHYD, basename
    
#==============================================================================
# DETERMINE CALIBRATION VALUE + ELECTRONICS CORRECTION
#==============================================================================

def CalibrationArray(file_electronics, file_ijkDOM, file_ijkREF, lowcut, 
                     highcut, cal_freq, NFFT, Fs_DOM, Fs_REF, cal_min, 
                     cal_max):

    print('\nGOING TO MAKE CALIBRATION ARRAY')

    #elec_freq, elec_ratio = Elec_ratio(file_electronics, cal_freq, NFFT, 
    #                                   lowcut, highcut)
#   dimension less, since it's a ratio

    ijkDOM, basenameDOM = get_cut_filt_IJkfile(file_ijkDOM, Fs_DOM, lowcut, 
                                               highcut, NFFT)
    ijkREF, basenameREF = get_cut_filt_IJkfile(file_ijkREF, Fs_REF, lowcut, 
                                               highcut, NFFT)
#==============================================================================
#==============================================================================
#==============================================================================
# # #                                                                                   <== hier spelen
#==============================================================================
#==============================================================================
#==============================================================================

    # PSD in dB re amp^2
    t1 = MakeTimeSpan(len(ijkDOM), int(25e6/128))
#    plt.plot((t1-0.01)/1e-3, ijkDOM, linewidth=1, label = 'DOM piezo')                                 # <==  TIME [MS]
    plt.plot((t1-0.01)*1500, ijkDOM, linewidth=1, label = 'DOM piezo')
    t2 = MakeTimeSpan(len(ijkREF), int(2e5))
#    plt.plot((t2-0.010015)/1e-3, ijkREF, color = 'r', linewidth=1, label = 'Reference hydrophone')     # <==  TIME [MS]
    plt.plot((t2-0.010015)*1500, ijkREF, color = 'r', linewidth=1, label = 'Reference hydrophone')
    plt.tick_params(axis='both', labelsize=15)
    plt.title('Reflection pattern wavelet 20kHz', fontsize=24)
    plt.ylabel("Amplitude [Pa]", fontsize=18)
    plt.legend(fontsize=18)
    plt.autoscale(tight = True)
    plt.grid()
    plt.grid('on', 'minor')
#    plt.xlabel('Time [ms]', fontsize=18)                                                                # <==  TIME [MS]
    plt.xlabel('Distance [m]', fontsize=18)
    plt.tight_layout()
    plt.show()
    plt.clf()

    raise sys.exit()


    def align_yaxis(ax1, v1, ax2, v2):
        """adjust ax2 ylimit so that v2 in ax2 is aligned to v1 in ax1"""
        _, y1 = ax1.transData.transform((0, v1))
        _, y2 = ax2.transData.transform((0, v2))
        inv = ax2.transData.inverted()
        _, dy = inv.transform((0, 0)) - inv.transform((0, y1-y2))
        miny, maxy = ax2.get_ylim()
        ax2.set_ylim(miny+dy, maxy+dy)

    fig, ax1 = plt.subplots()
    t1 = MakeTimeSpan(len(ijkDOM), int(25e6/128))
    ax1.plot(t1-0.01, ijkDOM, 'b-', linewidth=1)
    ax1.set_xlabel('time (s)', fontsize = 18)
    # Make the y-axis label, ticks and tick labels match the line color.
    ax1.set_ylabel('ADC', color='b')
    ax1.set_ylabel('amplitude [Pa]\nDOM piezo', color='b', fontsize = 18)
    ax1.tick_params('y', colors='b')
    ax1.tick_params(axis='both', labelsize=15)    
    ax1.grid()
    ax1.grid('on','both','minor')

    ax2 = ax1.twinx()
    t2 = MakeTimeSpan(len(ijkREF), int(2e5))
    ax2.plot(t2-0.010015, ijkREF, 'r-', linewidth=1)
    ax2.set_ylabel('amplitude [Pa]\nReference hydrophone ', color='r', fontsize = 18)
    ax2.tick_params('y', colors='r')
    ax2.tick_params(axis='both', labelsize=15)
    ax2.set_title('Reflection pattern wavelet 20kHz', fontsize=24)
    ax2.grid()
    ax2.grid('on','both', 'minor')
    
    #fig.suptitle('Reflection pattern wavelet 20kHz', fontsize=24)
    fig.tight_layout()
    align_yaxis(ax1, 0, ax2, 0)
    plt.show()
    plt.clf()

    raise sys.exit()    


    # verschil in sample frequentie leidt tot scheef lopen van de data
    # omdat data geknipt is op aantal sample punten.

#==============================================================================
#==============================================================================
#==============================================================================
# # # 
#==============================================================================
#==============================================================================
#==============================================================================
    DOM_val_X = IJk_value(ijkDOM, basenameDOM, int(25e6/128), NFFT, cal_min, 
                          cal_max)  # in [X]
    REF_val_Pa = IJk_value(ijkREF, basenameREF, int(2e5), NFFT, cal_min, 
                           cal_max) # in [Pa]

    sensitivity_value = (DOM_val_X / REF_val_Pa)   # in [X]/[Pa]
    dist_cor = distance_correction(config = 2)
    calibration_value = sensitivity_value / dist_cor

    print '- Used values used for calibration:'
    print '   Sensitivity value (DOM/REF) = %s' %sensitivity_value
    print '   Distance correction value   = %s' %dist_cor
    print '   Resulting calibration value = %s' %calibration_value

    return elec_freq, elec_ratio, calibration_value

#==============================================================================
# MAIN
#==============================================================================

def main():
    #start timer
    startTime = datetime.now()

    # input values for calubration:
# Specify type of calibration for savingname:
    cal_type        = 'wav10kHz'
# electronics file to remove 20dB offset:
    file_electronics= '20170419_6_white_noise_0.9Vrms.WAV'           # .wav
# file to determine calibration value Hydrophone:                                    <== hier bestandsnaam aangepast!
    file_ijkDOM     = 'c_wav10kHz_20170310_config2_DOM_8_v1_wavelet_20kHz_bp1-80.dat'  # .dat
# file to determine calibration value Reference:                                     <== hier bestandsnaam aangepast!
    file_ijkREF     = 'Time_Conf2Domwavelet20kHzv1.dat'              # .dat
# lower bound bandpass filter:
    lowcut          = 1000                                           # [Hz]
# upper bound bandpass filter:
    highcut         = 80000                                          # [Hz]
# frequency at which calibration is being done:
    cal_freq        = 10000                                          # [Hz]
# amount of samples to cut around pulse:
    samples_to_cut  = 1024                                           # [samp]
# samples to take (sub)PSD over:
    NFFT            = samples_to_cut # do NOT change this!           # [samp]
# sample frequency of DOM:
    Fs_DOM          = int(25e6/128)                                  # [samp/s]
# sample frequency of REF:
    Fs_REF          = int(2e5)                                       # [samp/s]
# lowerbound frequency for cummulative, to determine calibration value:
    cal_min         = 3500                                           # [Hz]
# upperbound frequency for cummulative, to determine calibration value:
    cal_max         = 22500                                          # [Hz]

    #get electronics ratio, and calibration value
    elec_freq, elec_ratio, calibration_value = CalibrationArray(
        file_electronics, file_ijkDOM, file_ijkREF, lowcut, highcut, cal_freq, 
        NFFT, Fs_DOM, Fs_REF, cal_min, cal_max)
    # due to different Fs, elec_freq is slightly different from freq --> 
    # results in about 1.6% deviation in samples



    #end timer
    print datetime.now() - startTime

if __name__ == "__main__":
    sys.exit(main())
