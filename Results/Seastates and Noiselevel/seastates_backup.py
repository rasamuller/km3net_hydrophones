# -*- coding: utf-8 -*-
"""
Created on Fri Jun 09 01:40:27 2017

@author: mullerrs
"""

import sys
import numpy as np
import struct
import os
import math
import scipy.io
import matplotlib.pyplot as plt
from datetime import datetime
from scipy.signal import butter, lfilter
from numpy.linalg import norm as Norm
from numpy.fft import fft


def MakeFreqSpan(Fs, NFFT):
    return np.linspace(0, NFFT/2-1, NFFT/2)*float(Fs)/float(NFFT)

def MakeSpectrum(Fs, NFFT, data, basename):

    length     = len(data)
    PSD        = np.zeros((NFFT/2, 1), dtype=float)
    freq       = MakeFreqSpan(Fs, NFFT)
    Segment    = int(length/NFFT)
    if 'sweep' in basename:
        wnd    = np.ones(NFFT)
        print 'rectangular window used'
    else:
        wnd    = np.hanning(NFFT)
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

def main():
    Fs       = int(25e6/128)
    NFFT     = 1024
    freq     = np.linspace(0, NFFT/2-1, NFFT/2)*float(Fs)/float(NFFT)
    ss       = [1., 2., 3., 4., 5., 6., 7., 8., 9., 10.]
    noise    = np.genfromtxt('c_wav10kHz_20170309_config2_DOM_15_v2_ruis_bp1-80.dat')
    PSDnoise = MakeSpectrum(Fs, NFFT, noise, 'c_wav10kHz_20170309_config2_DOM_15_v2_ruis_bp1-80.dat')[0]

    # KUNDSEN - Kundsen et al., 1948
    # ss0 impossible since log(0) = - inf
    plt.plot(freq, 20*np.log10(np.sqrt(PSDnoise)/(1e-6))-10*np.log10(Fs/float(NFFT)), linewidth = 1.5, label = 'noise level DOM')
    plt.tick_params(axis='both', labelsize=15)
    plt.autoscale(tight = True)
    plt.grid()
    plt.grid('on', 'minor')
    plt.tight_layout()

    for i in ss:
        plt.plot(freq, 56. + 19.* np.log10(i) - 17. * np.log10(freq/1000.), '--', label='seastate %s' %int(i))#str(i))
        plt.legend(bbox_to_anchor=(1.04, 1.02), fontsize = 'small')
    plt.title("\n\n Noise DOM piezo v.s.\nambient noise level per seastate", fontsize = 24)
    plt.suptitle('\n based on Kundsen et al., 1948', fontsize = 9)
    plt.xlabel("Frequency (Hz)", fontsize = 18)
    plt.ylabel("dB re 1 $\mu$Pa/$\sqrt{Hz}$", fontsize = 18)
    #plt.xlim([1000,25000])
    plt.ylim([10, 80])
    plt.tight_layout()
    #plt.text(1.3, 1.8,'NL = 56 + 19 log (ss) - 17 log (f)', horizontalalignment='center', verticalalignment='center', fontsize = 'medium', family = 'serif')
#    plt.savefig("seastates_knudsen_1984.png", dpi=200)
    plt.show()
    plt.clf()

if __name__ == "__main__":
    sys.exit(main())