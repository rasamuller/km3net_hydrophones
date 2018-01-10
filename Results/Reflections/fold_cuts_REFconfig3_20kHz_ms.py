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
    time = np.linspace(0, N-1, N)/Fs
    for segment in range(int(len(cuttedfilename)/N)):
        segment = segment + 1
        if len(cuttedfilename[segment*N: (segment+1)*N]) == len(time):
            plt.plot(time/(1e-3)-57.5, cuttedfilename[segment*N: (segment+1)*N]/(1e5))
        break
    plt.title('\n\nSingle pulse of Reference hydrophone', fontsize=24)
    plt.suptitle('used: %s' %(basename), fontsize=9)
    plt.tick_params(axis='both', labelsize=15)
    plt.xlabel('Time [ms]', fontsize=18)
    #plt.xlabel('Distance [m]', fontsize=18)
    plt.autoscale(tight = True)
    plt.grid()
    plt.grid(True, zorder = 1)
    #plt.grid('on', 'minor')
#    plt.legend()
    plt.tight_layout()
    #plt.savefig(basename +'.png')
    plt.show()
    plt.clf()

def main():
    N = 50000#4096
    Fs = 2e5

    data, basename = Import_data('Time_Conf3Barewavelet20kHzv1.dat')
    data, basename = butter_bandpass_filter(data, basename, 1000, 80000, Fs, order=5)

    PlotFoldedCuts(data, basename, N, Fs)


if __name__ == "__main__":
    sys.exit(main())
