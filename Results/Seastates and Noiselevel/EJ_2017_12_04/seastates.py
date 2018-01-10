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

def main():

    filename = 'c_wav10kHz_20170309_config2_DOM_15_v2_ruis_bp1-80_wnd_hann_spectrumPSD.dat'
#    filename = 'c_wav10kHz_20170309_config2_DOM_15_v2_ruis_bp1-80_wnd_hann_spectrumPSD_ELcor.dat'

    NFFT     = 1024
    Fs       = 25e6/128
    freq     = np.linspace(0, NFFT/2-1, NFFT/2)*float(Fs)/float(NFFT)
    PSD      = np.genfromtxt(filename)
    basename = os.path.splitext(os.path.basename(filename))[0] 
    #ss       = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    ss       = [6, 4, 3, 2, 1, 0.5]
    ws       = ['28-33', '17-21', '11-16', '7-10', '4-6', '1-3']

    
    # NOISE LEVEL DOM
    plt.plot(freq, 20*np.log10(np.sqrt(PSD)/(1e-6))-10*np.log10(Fs/float(NFFT)), linewidth = 1.5, label = 'Noise level \nDOM piezo')

    # KUNDSEN - Kundsen et al., 1948
    # ss0 impossible since log(0) = - inf
    for i in ss:
        freq1 = freq[:262]
        plt.plot(freq1, 56. + 19.* np.log10(i) - 17. * np.log10(freq1/1000.), '-', label='seastate %s' %i)#str(i))
        freq2 = freq[262:]
        plt.plot(freq2, 56. + 19.* np.log10(i) - 17. * np.log10(freq2/1000.), '--k')
        plt.legend(bbox_to_anchor=(1.04, 1.02), fontsize = 'small')

    plt.xlim(0,262)
    plt.title("\n\n Noise DOM piezo v.s.\nambient noise level per seastate", fontsize = 24)
    plt.suptitle('\n used: %s\n based on Kundsen et al., 1948' %(basename), fontsize = 9)
    plt.xlabel("Frequency (Hz)", fontsize = 18)
    plt.ylabel("dB re 1 $\mu$Pa/$\sqrt{Hz}$", fontsize = 18)
    plt.tick_params(axis='both', labelsize=15)
    plt.autoscale(tight = True)
    plt.grid()
    plt.grid('on', 'minor')
   # plt.ylim([10, 80])
    plt.tight_layout()
#    plt.savefig("seastates_knudsen_1984.png", dpi=200)
    plt.show()
    plt.clf()

if __name__ == "__main__":
    sys.exit(main())