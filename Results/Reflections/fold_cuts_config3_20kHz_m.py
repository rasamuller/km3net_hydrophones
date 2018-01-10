# -*- coding: utf-8 -*-
"""
Created on Fri May 12 12:04:43 2017

@author: mullerrs
"""

import sys
import numpy as np
import struct
import os
import math
import scipy.io
from SignalProc_BAD import *
import matplotlib.pyplot as plt
from datetime import datetime
from scipy.signal import butter, lfilter


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

def Import_data(filename):
    filename = '.\\' + filename
    print filename
    basename = os.path.splitext(os.path.basename(filename))[0]
    return np.genfromtxt(filename), basename

def PlotFoldedCuts(cuttedfilename, basename, N, Fs):
    basename = basename + '_fold'
    time = np.linspace(0, N-1, N)/Fs*1500
    for segment in range(int(len(cuttedfilename)/N)):
        segment = segment + 1
        if len(cuttedfilename[segment*N: (segment+1)*N]) == len(time):
            # 4 wavelets naast elkaar            
            #plt.plot(time[450:-150], cuttedfilename[segment*N+450: (segment+1)*N-150]/(1e5))
            # 1 wavelet wat langer
            plt.plot(time-84.2, cuttedfilename[segment*N: (segment+1)*N]/(1e5))
        break
    plt.title('\n\nSingle pulse of bare piezo \nin waterproof box', fontsize=24)
    plt.suptitle('used: %s' %(basename), fontsize=9)
    plt.tick_params(axis='both', labelsize=15)
    #plt.xlabel('Time [ms]', fontsize=18)
    plt.xlabel('Distance [m]', fontsize=18)
    plt.autoscale(tight = True)
    plt.grid()
    plt.grid(True, zorder = 1)
    #plt.grid('on', 'minor')
#    plt.legend()
    plt.tight_layout()
    #plt.savefig(basename +'.png')
    plt.show()
    plt.clf()

def Array_mean(cuttedfilename, basename, N):
    basename = basename + '_fold_mean'
    ADD = np.zeros(N)
    for segment in range(int(len(cuttedfilename)/N)-1):
        segment = segment + 1
        ADD = ADD + cuttedfilename[segment*N: (segment+1)*N]
        plt.plot(cuttedfilename[segment*N: (segment+1)*N], 'g.')
    MEAN = ADD / (int(len(cuttedfilename)/N))
    #np.savetxt(basename + '.dat', MEAN)
    plt.plot(MEAN, 'r-', label = 'mean')
    plt.xlabel('samples -->')
    plt.ylabel('ADC count --> ')
    plt.title(basename)
    plt.legend()
    #plt.savefig(basename +'.png')
    plt.show()
    plt.clf()
    return MEAN

def main():
    N = 50000#4096
    Fs = 25e6/128
#    Fs = 2e5

    #data, basename = Import_data('20170310_config2_DOM_6_v1_wavelet_10kHz_bpfilter_beforecutting__pulses_cut.dat')
    #data, basename = Import_data('20170310_config2_DOM_8_v1_wavelet_20kHz_bpfilter_beforecutting__pulses_cut.dat')
#    data, basename = Import_data('20170310_config2_DOM_6_v1_wavelet_10kHz.dat')    
    data, basename = Import_data('20170310_config3_Bare_18_v1_wavelet_20kHz.dat')    
    #data, basename = Import_data('Time_Conf2Domwavelet10kHzv1_bpfilter_beforecutting__pulses_cut.dat')
    #data, basename = Import_data('Time_Conf2Domwavelet20kHzv1_bpfilter_beforecutting__pulses_cut.dat')

    data, basename = butter_bandpass_filter(data, basename, 1000, 80000, Fs, order=5)

    PlotFoldedCuts(data, basename, N, Fs)
#    mean = Array_mean(data, basename, N)
#    plt.plot(mean)
#    plt.show()
#    plt.clf()
#    
#    ratio, ratio_basename = Import_data('20170419_6_white_noise_0.9Vrms_RATIO_nomalized.dat')    
#    
#    f, PSD_mean = signal.welch(mean, int(Fs), nperseg=1024)
#    f, PSD_dat = signal.welch(data, int(Fs), nperseg=1024)
#    plt.plot(f, 20*np.log10(np.sqrt(PSD_mean)), 'r.', label = 'mean')
#    plt.plot(f, 20*np.log10(np.sqrt(PSD_dat)), 'g.', label = 'dat')
#    plt.plot(f, 20*np.log10(np.sqrt(PSD_mean/ratio)), 'r-', label = 'mean/ratio')
#    plt.plot(f, 20*np.log10(np.sqrt(PSD_dat/ratio)), 'g-', label = 'dat/ratio')
#    plt.legend()
#    plt.show()
#    plt.clf()


#==============================================================================
# handmatig --> MAG WEG
#==============================================================================
#    data10 = np.genfromtxt('20170310_config2_DOM_6_v1_wavelet_10kHz_pulses_cut.dat')
#    for segment in range(int(len(data10)/N)):
#        segment = segment + 1
#        plt.plot(data10[segment*N : (segment+1)*N])
#        plt.title('DOM 10kHz, not filtered')
#    plt.show()
#    plt.clf()
#
#    fdata10 = np.genfromtxt('20170310_config2_DOM_6_v1_wavelet_10kHz_bpfilter_beforecutting__pulses_cut.dat')
#    for segment in range(int(len(fdata10)/N)):
#        segment = segment + 1
#        plt.plot(fdata10[segment*N : (segment+1)*N])
#        plt.title('DOM 10kHz, filtered')
#    plt.show()
#    plt.clf()
#
#    data20 = np.genfromtxt('20170310_config2_DOM_8_v1_wavelet_20kHz_pulses_cut.dat')
#    for segment in range(int(len(data20)/N)):
#        segment = segment + 1
#        plt.plot(data20[segment*N : (segment+1)*N])
#        plt.title('DOM 20kHZ, not filtered')
#    plt.show()
#    plt.clf()
#
#    fdata20 = np.genfromtxt('20170310_config2_DOM_8_v1_wavelet_20kHz_bpfilter_beforecutting__pulses_cut.dat')
#    for segment in range(int(len(fdata20)/N)):
#        segment = segment + 1
#        plt.plot(fdata20[segment*N : (segment+1)*N])
#        plt.title('DOM 20kHz, filtered')
#    plt.show()
#    plt.clf()
#
#    ref10 = np.genfromtxt('Time_Conf2Domwavelet10kHzv1_pulses_cut.dat')
#    for segment in range(int(len(ref10)/N)):
#        segment = segment + 1
#        plt.plot(ref10[segment*N : (segment+1)*N])
#        plt.title('REF 10kHz, not filtered')
#    plt.show()
#    plt.clf()
#
#    fref10 = np.genfromtxt('Time_Conf2Domwavelet10kHzv1_bpfilter_beforecutting__pulses_cut.dat')
#    for segment in range(int(len(fref10)/N)):
#        segment = segment + 1
#        plt.plot(fref10[segment*N : (segment+1)*N])
#        plt.title('REF 10kHz, filtered')
#    plt.show()
#    plt.clf()
#
#    ref20 = np.genfromtxt('Time_Conf2Domwavelet20kHzv1_pulses_cut.dat')
#    for segment in range(int(len(ref20)/N)):
#        segment = segment + 1
#        plt.plot(ref20[segment*N : (segment+1)*N])
#        plt.title('REF 20kHZ, not filtered')
#    plt.show()
#    plt.clf()
#
#    fref20 = np.genfromtxt('Time_Conf2Domwavelet20kHzv1_bpfilter_beforecutting__pulses_cut.dat')
#    for segment in range(int(len(fref20)/N)):
#        segment = segment + 1
#        plt.plot(fref20[segment*N : (segment+1)*N])
#        plt.title('REF 20kHz, filtered')
#    plt.show()
#    plt.clf()

if __name__ == "__main__":
    sys.exit(main())
