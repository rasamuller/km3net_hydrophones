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
def MakeSpectrum(Fs, NFFT, data, basename):

    length     = len(data)
    PSD        = np.zeros((NFFT/2, 1), dtype=float)
    freq       = np.linspace(0, NFFT/2-1, NFFT/2)*float(Fs)/float(NFFT)
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

for fn in os.listdir('.'):
#    if "spectrumPSD.dat" in fn:
    l = ("20170309_config2_DOM_21_sweep10-60kHz_1Vpp.dat", "c_wav10kHz_20170309_config2_DOM_21_sweep10-60kHz_1Vpp_bp1-80_wnd_rect_spectrumPSD.dat", "c_wav10kHz_20170309_config2_DOM_21_sweep10-60kHz_1Vpp_bp1-80_wnd_rect_spectrumPSD_ELcor.dat")
    for i in l:

        if i in fn:
            Fs          = 25e6/128
            NFFT        = 1024
            freq        = np.linspace(0, NFFT/2-1, NFFT/2)*float(Fs)/float(NFFT)

            if 'Vpp.dat' in i:
                data = np.genfromtxt(fn)
                basename    = os.path.splitext(os.path.basename(fn))[0]
                PSD, basename = MakeSpectrum(Fs, NFFT, data, basename)
                plt.plot(freq, (10*np.log10(abs(PSD)/(1e-6)**2))-(10*np.log10((25e6/128)/1024))-60, label = "not calib, \n but taken down")

            else:
                PSD        = np.genfromtxt(fn)
                basename    = os.path.splitext(os.path.basename(fn))[0]
                if 'ELcor.dat' in i:
                    labl = "calib, \n  Elec corr"
                elif 'PSD.dat' in i:
                    labl = "calib, \n not Elec corr"
                #labl        = basename.split('_')[3] # Voor DOM: welke config
                #typ         = basename.split('_')[7]
               
                plt.plot(freq, (10*np.log10(abs(PSD)/(1e-6)**2))-(10*np.log10((25e6/128)/1024)), label = '%s' %labl, linewidth=1) # dB re $\mu$Pa
            plt.tick_params(axis='both', labelsize=15)
            plt.title('Sweep 10 - 60 kHz, config 2, 1Vpp', fontsize=24)
            plt.legend(fancybox = True, shadow = True, loc = 'best')
            plt.ylabel('Amplitude [dB re $\mu$Pa/$\sqrt{Hz}$]', fontsize=18)
            plt.xlabel('Freqency [Hz]', fontsize=18)
            plt.savefig('noise_spectra_psD_corr')
plt.grid()
plt.grid('on', axis = 'minor')
plt.show()
plt.clf()