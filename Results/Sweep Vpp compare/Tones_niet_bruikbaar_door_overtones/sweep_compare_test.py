# -*- coding: utf-8 -*-
"""
Created on Thu May 25 22:25:27 2017

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

Fs = int(25e6/128)
Fs = int(2e5)
audioFs = 15.25
show = False
save_png = False
save_txt = False


#==============================================================================
# soort van de main
#==============================================================================
for fn in os.listdir('.'):
    if "ELcor.dat" and "_10kHz_" in fn:
#    if "PSD.dat" in fn:
    #l = ("c_wav10kHz_20170309_config2_DOM_21_sweep10-60kHz_1Vpp_bp1-80_wnd_rect_spectrumPSD.dat", "c_wav10kHz_20170309_config2_DOM_21_sweep10-60kHz_1Vpp_bp1-80_wnd_rect_spectrumPSD_ELcor.dat")
    #for i in l:
    #    if i in fn:
        PSD        = np.genfromtxt(fn)
        basename    = os.path.splitext(os.path.basename(fn))[0]
        labl        = basename.split('_')[3] # Voor DOM: welke config
        typ         = basename.split('_')[8]
        dF          = int(25e6/128)/1024
    
            # create black dots in the graph to show which poins are used
            # for the linearity plot.
        plt.plot(10000 ,10*np.log10(abs(PSD[10000/dF])/(1e-6)**2)-(10*np.log10(dF)), 'ok')
        plt.plot(20000 ,10*np.log10(abs(PSD[20000/dF])/(1e-6)**2)-(10*np.log10(dF)), 'ok')
        plt.plot(30000 ,10*np.log10(abs(PSD[30000/dF])/(1e-6)**2)-(10*np.log10(dF)), 'ok')
        plt.plot(40000 ,10*np.log10(abs(PSD[40000/dF])/(1e-6)**2)-(10*np.log10(dF)), 'ok')
        plt.plot(50000 ,10*np.log10(abs(PSD[50000/dF])/(1e-6)**2)-(10*np.log10(dF)), 'ok')
        plt.plot(60000 ,10*np.log10(abs(PSD[60000/dF])/(1e-6)**2)-(10*np.log10(dF)), 'ok')
        plt.plot(70000 ,10*np.log10(abs(PSD[70000/dF])/(1e-6)**2)-(10*np.log10(dF)), 'ok')

        Fs          = 25e6/128
        NFFT        = 1024
        freq        = np.linspace(0, NFFT/2-1, NFFT/2)*float(Fs)/float(NFFT)

        plt.plot(freq, (10*np.log10(abs(PSD)/(1e-6)**2))-(10*np.log10((25e6/128)/1024)), label = '%s' %typ, linewidth=1) # dB re $\mu$Pa
        plt.tick_params(axis='both', labelsize=15)
        plt.title('Sweep 10-60 kHz', fontsize=24)
        plt.legend(fancybox = True, shadow = True, loc = 'best')
        plt.ylabel('Amplitude [dB re $\mu$Pa/$\sqrt{Hz}$]', fontsize=18)
        plt.xlabel('Freqency [Hz]', fontsize=18)
        plt.xlim(10000, 75000)
        #plt.savefig('sweep_spectra_psD_corr_noelcor')
plt.grid()
plt.grid('on', axis = 'minor')
plt.show()
plt.clf()