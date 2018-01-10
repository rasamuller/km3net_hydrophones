# -*- coding: utf-8 -*-
"""
Created on Thu Jun 15 14:36:18 2017

@author: mullerrs
"""
import numpy as np
import sys
import os
from matplotlib import pyplot as plt

freq = (10000., 20000., 30000., 50000., 60000.)     # in Hz
freq = np.array(freq)
TVR_dB = (110., 124., 128., 135., 141.) # in dB re micro Pa/V
TVR_dB = np.array(TVR_dB)

TVR_muPaV = 10**(TVR_dB/20) # in micro Pa/V

Vpp = [1, 3, 5, 7, 10]

dist_corr = (4*np.pi*(0.97**2))

print '======================'
print 'distance correction = ', dist_corr
print '======================'
print 'For: x Vpp, transmitter produces: x [Pa] for:'
print '          10kHz         20kHz        30kHz        50kHz        60kHz     '

for i in Vpp:
    TVR_muPa = i*TVR_muPaV
    TVR_Pa = TVR_muPa/1e6
    print '%sVpp'%i, TVR_Pa

print '\nFor: x Vpp, DOM should measure: x [Pa] for:'
print '          10kHz         20kHz        30kHz        50kHz        60kHz     '

for i in Vpp:
    TVR_muPa = i*TVR_muPaV
    TVR_Pa = TVR_muPa/1e6
    print '%sVpp'%i, TVR_Pa/dist_corr


dist_corr_REF = (4*np.pi*(0.45**2))

#################
print '\nFor: x Vpp, REF should measure: x [Pa] for:'
print '          10kHz         20kHz        30kHz        50kHz        60kHz     '

for i in Vpp:
    TVR_muPa = i*TVR_muPaV
    TVR_Pa = TVR_muPa/1e6
    print '%sVpp'%i, TVR_Pa/dist_corr_REF
################


#for fn in os.listdir('.'):
#    if '_50kHz' in fn and '.dat' in fn:
#        f = 50000
#        maxrows = 10000
#        non_calibrated_data = np.genfromtxt(fn, usecols=0, max_rows = maxrows)
#
#        basename = os.path.splitext(os.path.basename(fn))[0]
#        print basename
#        vpp = basename.split('_')[-1]
#
#        Fs = int(25e6/128)
#        length = len(non_calibrated_data)
#
#        # mean max
#        ma = np.max(non_calibrated_data)
#        ma = ma - 0.03*ma
#        # mean min
#        mi = np.min(non_calibrated_data)
#        mi = mi - 0.03*mi
#        #amp
#        AMPpp = ma - mi
#
#        time = np.linspace(0, length-1, length)/Fs
#
#        plt.plot(time, non_calibrated_data)
#        plt.plot(time, np.ones(len(time))*ma, 'r--')
#        plt.plot(time, np.ones(len(time))*mi, 'r--')
#        plt.title('%skHz tone, %s \nAMPpp = %s' %(int(f/1000), vpp, AMPpp))
#        plt.xlabel('Time (s)')
#        plt.ylabel('Amplitude')
#        plt.tight_layout()
#        plt.grid()
#        plt.grid('on', 'minor')
#        plt.savefig(basename+'_data%s.png' %maxrows)
#        plt.show()
#        plt.clf()

#==============================================================================
# 10 kHZ
#==============================================================================
Pa10  = [0.027, 0.080, 0.134, 0.187, 0.267]
Pa10  = np.array(Pa10)
AMP10 = [54026.79, 163370.31, 263161.97, 369736.84, 527364.75]
AMP10 = np.array(AMP10)

#==============================================================================
# 30 kHz
#==============================================================================

Pa30  = [0.212, 0.637, 1.062, 1.487, 2.124]
Pa30  = np.array(Pa30)
AMP30 = [555855.59, 1627306.92, 2677850.87, 3723636.0, 5274860.97]
AMP30 = np.array(AMP30)

#==============================================================================
# 50kHz
#==============================================================================

Pa50  = [0.476, 1.427, 2.378, 3.329, 4.756]
Pa50  = np.array(Pa50)
AMP50 = [1812232.57, 5464798.61, 8634992.38, 12569091.22, 16273898.0]
AMP50 = np.array(AMP50)

#==============================================================================
# PLOT
#==============================================================================

plt.plot(Vpp, AMP10/Pa10, 'o-', linewidth = 1, label = '10kHz')
plt.plot(Vpp, AMP30/Pa30, 'o-', linewidth = 1, label = '30kHz')
plt.plot(Vpp, AMP50/Pa50, 'o-', linewidth = 1, label = '50kHz')
plt.title('DOM Amp/Pa per transmtter Vpp', fontsize = 24)
plt.tick_params(axis='both', labelsize=15)
plt.xlabel('Vpp', fontsize = 18)
plt.ylabel('AMP/Pa', fontsize = 18)
plt.legend()
plt.tight_layout()
plt.savefig('DOM_tones_compTransm')
plt.show()
plt.clf()

