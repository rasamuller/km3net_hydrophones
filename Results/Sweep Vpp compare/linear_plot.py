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
    if "ELcor.dat" in fn:
    #if "PSD.dat" in fn:
    #l = ("c_wav10kHz_20170309_config2_DOM_21_sweep10-60kHz_1Vpp_bp1-80_wnd_rect_spectrumPSD.dat", "c_wav10kHz_20170309_config2_DOM_21_sweep10-60kHz_1Vpp_bp1-80_wnd_rect_spectrumPSD_ELcor.dat")
    #for i in l:
    #    if i in fn:
        PSD        = np.genfromtxt(fn)
        basename    = os.path.splitext(os.path.basename(fn))[0]
        labl        = basename.split('_')[3] # Voor DOM: welke config
        typ         = basename.split('_')[7]
        dF          = int(25e6/128)/1024
        Fs          = 25e6/128
        NFFT        = 1024
        freq        = np.linspace(0, NFFT/2-1, NFFT/2)*float(Fs)/float(NFFT)

        f   = [10000, 20000, 40000, 60000]
        lin = [PSD[10000/dF], PSD[20000/dF], PSD[40000/dF], PSD[60000/dF]]
        lin = np.array(lin)

        if typ == '1Vpp':
            Vpp1 = lin
            print typ, lin
        elif typ == '3Vpp':
            Vpp3 = lin
            print typ, lin
        elif typ == '5Vpp':
            Vpp5 = lin
            print typ, lin
        elif typ == '7Vpp':
            Vpp7 = lin
            print typ, lin
        elif typ == '10Vpp':
            Vpp10 = lin
            print typ, lin

for i in range(len(Vpp10)):
    print f[i]
    print Vpp1[i]
    print Vpp3[i]
    print Vpp5[i]
    print Vpp7[i]
    print Vpp10[i]
    vpp_y = [Vpp1[i], Vpp3[i], Vpp5[i], Vpp7[i], Vpp10[i]]
    vpp_x = [1, 3, 5, 7, 10]

    # determine parameters for linear fit: y=ax+b
    a = [(Vpp10[i]-Vpp1[i])/(10-1)]
    a = np.array(a)
    b = Vpp10[i]-a*10
    print '----b----:', b
    b = np.array(b)

    labl = f[i]/1000

    plt.plot(vpp_x, vpp_y, 'o', label = '%s kHz' %labl, linewidth=1)
    plt.plot(vpp_x, vpp_y, '--k', linewidth=1)

## LINEAR FIT
#    plt.plot(vpp_x, a*vpp_x+b, '--r', linewidth = 1)
## just repeat the last one for the lable to appear only once
#plt.plot(vpp_x, a*vpp_x+b, '--r', linewidth = 1, label = 'linear fit')

## y=ax^k FIT
#for k in range(1):
#    k = k/10. + 1.5
#    print 'k = ', k
#    k = 2
#    for a in range(10):
#        a = a/10.
#        #a=0.1
#        print 'a = ', a
#        X = []
#        Y = []
#        for x in range(10):
#            X.append(x)
#            Y.append(a*x**k)
#        plt.plot(X, Y, label = a)


#plt.yscale("log")
#plt.xscale("log")
plt.tick_params(axis='both', labelsize=15)
plt.title('Sweeps 10-60 kHz', fontsize=24)
plt.legend(fancybox = True, shadow = True, loc = 'best')
plt.ylabel('Amplitude [Pa]', fontsize=18)
plt.xlabel('Input voltage source [Vpp]', fontsize=18)
plt.grid()
plt.grid('on', axis = 'minor')
#plt.savefig('linear_behaviour_Vpp')
plt.show()
plt.clf()