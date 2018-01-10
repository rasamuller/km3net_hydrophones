# -*- coding: utf-8 -*-
"""
Created on Wed Apr 12 15:17:57 2017

@author: mullerrs
"""

#import pywt
import matplotlib.pyplot as plt
import numpy as np
import math

#
####################
####################
####################
##to the axes bounding box, with 0,0 being the lower left of the axes and 1,1 the upper right.
#import matplotlib.pyplot as plt
#import matplotlib.patches as patches
#
## build a rectangle in axes coords
#left, width = .25, .5
#bottom, height = .25, .5
#right = left + width
#top = bottom + height
#
#fig = plt.figure()
#ax = fig.add_axes([0,0,1,1])
#
## axes coordinates are 0,0 is bottom left and 1,1 is upper right
#p = patches.Rectangle(
#    (left, bottom), width, height,
#    fill=False, transform=ax.transAxes, clip_on=False
#    )
#
#ax.add_patch(p)
####################
####################
####################

#    plt.plot(freq, 10.*np.log((freq/1000.)**(-5/3)) + 94.5 + 30*np.log(ss[i-1] + 1), color = 'black')
#    plt.plot(freq, np.sqrt(10.*np.log(freq**(-5/3)) + 94.5 + 30*np.log(ns[i-1] + 1)))


Fs       = 25e6/128
NFFT     = int(Fs/10.)
freq     = np.linspace(0, NFFT/2-1, NFFT/2)*float(Fs)/float(NFFT)
ss       = [1., 2., 3., 4., 5., 6.]


######################################
######################################
#####                            #####
##### - Hildebrand, J.A., 2009 - #####
#####                            #####
######################################
######################################


################################################################################
#
## KUNDSEN - Kundsen et al., 1948
## ss0 impossible since log(0) = - inf
#for i in ss:
#    plt.plot(freq, 56. + 19.* np.log10(i) - 17. * np.log10(freq/1000.), label=str(i))
#    plt.legend(title='seastate', bbox_to_anchor=(1.04, 1.05), fontsize = 'small')
#    plt.title("AMBIENT NOISE LEVEL PER SEASTATE \n based on Kundsen et al., 1948")
#    plt.xlabel("Frequency (Hz)")
#    plt.ylabel("dB re 1 $\mu$Pa/$\sqrt{Hz}$")
#    plt.xlim([1000,25000])
#    plt.ylim([30, 70])
#    plt.tight_layout()
##plt.text(1.3, 1.8,'NL = 56 + 19 log (ss) - 17 log (f)', horizontalalignment='center', verticalalignment='center', fontsize = 'medium', family = 'serif')
#plt.savefig("seastates_knudsen_1984.png", dpi=200)
#plt.show()
#plt.clf()
#
################################################################################
#
# THERMAL NOISE - Mellen et al., 1952
plt.plot(freq, -15. + 20.* np.log10(freq/1000.))
for i in range(len(freq)):
    if freq[i]>50000:
        plt.plot(freq[i:len(freq)], -15. +20. * np.log10(freq[i:(len(freq))]/1000.), linewidth = 2, color = 'blue')
        break
plt.title("BACKGROUND LEVEL DUE TO THERMAL NOISE \n based on Mellen et al., 1952")
plt.xlabel("Frequency (Hz)")
plt.ylabel("dB re 1 $\mu$Pa")
plt.xlim([1000, 100000])
plt.ylim([0, 35])
#plt.text(1.3, 1.2,'NL = - 15 + 20 log (f)', horizontalalignment='center', verticalalignment='center', fontsize = 'medium', family = 'serif')
plt.tight_layout()
plt.savefig("background_level_due_to_thermal_noise.png", dpi = 200)
plt.show()
plt.clf()
#
################################################################################
#
## ATTENUATION - Ainslie & McColm, 1998
#S  = 0. #salinity %. 
#pH = 7. # <-- need to check
#T  = 18. # in #C
#z  = 0.002 # depth in km
#f1 = 0.78*np.sqrt(S/35)*math.exp(T/26)
#f2 = 42.*math.exp(T/17.)
## attenuation (dB/km), f in kHz
#A = 0.106*(f1*(freq/1000.)**2)/((freq/1000.)**2+f1**2)*math.exp((pH-8.)/0.56) + 0.52*(1+T/43)*(S/35)*(f2*(freq/1000.)**2)/((freq/1000.)**2+f2**2)*math.exp(-z/6) + 0.00049*(freq/1000.)**2*math.exp(-((T/27)+(z/17)))
#
#plt.plot(freq, A)
#plt.title("ATTENUATION COEFFICIENT (approximated) \n based on Ainslie & McColm, 1998")
#plt.xlabel("Frequency (Hz)")
#plt.ylabel("dB km$^{-1}$")
#plt.tight_layout()
#plt.text(0.5, 1.75,'Bath properties:\n S = %s '% S + '\n pH = %s'% pH + '\n T = %s '% T + '\n z = %s'% z, transform=ax.transAxes, horizontalalignment='left', verticalalignment='center', fontsize = 'small', family = 'serif', bbox = dict(facecolor='white', edgecolor='black', pad=10.0))
##plt.text(50.3, 2.2, r'$0.106 \frac{f_{1} f^{2}}{f^{2} + f_{1}^{2}} e^{(pH-8)/0.56} + 0.52(1+\frac{T}{43})(\frac{S}{35})\frac{f_{2} f^{2}}{f^{2}+f_{2}^{2}} e^{-z/6}$', horizontalalignment='center', verticalalignment='center', fontsize = 'large', family = 'serif')
##plt.text(1.3, 50.9, r'$+ 0.00049 f^{2} e^{-(T/27 + z/17)}$', horizontalalignment='center', verticalalignment='center', fontsize = 'large', family = 'serif')
#plt.savefig("attenuation_bath.png", dpi = 200)
#plt.show()
#plt.clf()
#
################################################################################
#
## ATTENUATION SALINITY- Ainslie & McColm, 1998
#S  = [0., 20., 40., 60., 80., 100.] #salinity %. 
#pH = 7. 
#T  = 18. # in #C
#z  = 0.002 # depth in km
#f2 = 42.*math.exp(T/17.)
## attenuation (dB/km), f in kHz
#
#for i in S:
#    f1 = 0.78*np.sqrt(i/35.)*math.exp(T/26.)
#    lw = 1
#    ls = '--'
#    if i == 0:
#        lw = 2
#        ls = '-'
#    plt.plot(freq, 0.106*(f1*(freq/1000.)**2)/((freq/1000.)**2+f1**2)*math.exp((pH-8.)/0.56) + 0.52*(1+T/43)*(i/35)*(f2*(freq/1000.)**2)/((freq/1000.)**2+f2**2)*math.exp(-z/6) + 0.00049*(freq/1000.)**2*math.exp(-((T/27)+(z/17))), label = str(i), linewidth = lw, linestyle = ls)
#    plt.title("ATTENUATION COEFFICIENT (approximated) \n based on Ainslie & McColm, 1998")
#    plt.legend(title='Salinity', bbox_to_anchor=(0.2, 0.95), fontsize = 'small')
#    plt.xlabel("Frequency (Hz)")
#    plt.ylabel("dB km$^{-1}$")
#    #plt.xlim([1000,25000])
#    #plt.ylim([30, 70])
#    plt.tight_layout()
##plt.text(50.3, 60.2, r'$0.106 \frac{f_{1} f^{2}}{f^{2} + f_{1}^{2}} e^{(pH-8)/0.56} + 0.52(1+\frac{T}{43})(\frac{S}{35})\frac{f_{2} f^{2}}{f^{2}+f_{2}^{2}} e^{-z/6}$', horizontalalignment='center', verticalalignment='center', fontsize = 'large', family = 'serif')
##plt.text(50.3, 50.9, r'$+ 0.00049 f^{2} e^{-(T/27 + z/17)}$', horizontalalignment='center', verticalalignment='center', fontsize = 'large', family = 'serif')
#plt.savefig("attenuation_bath_salinity.png", dpi = 200)
#plt.show()
#plt.clf()
#
################################################################################
#
## ATTENUATION pH- Ainslie & McColm, 1998
#S  = 0. #salinity %. 
#pH = [0., 3.5, 7., 10.5, 14]
#T  = 18. # in #C
#z  = 0.002 # depth in km
#f2 = 42.*math.exp(T/17.)
#f1 = 0.78*np.sqrt(S/35.)*math.exp(T/26.)
## attenuation (dB/km), f in kHz
#
#for i in pH:
#    lw = 1
#    ls = '--'
#    if i == 0:
#        lw = 2
#        ls = '-'
#    plt.plot(freq, 0.106*(f1*(freq/1000.)**2)/((freq/1000.)**2+f1**2)*math.exp((i-8.)/0.56) + 0.52*(1+T/43)*(S/35)*(f2*(freq/1000.)**2)/((freq/1000.)**2+f2**2)*math.exp(-z/6) + 0.00049*(freq/1000.)**2*math.exp(-((T/27)+(z/17)))
#, label = str(i), linewidth = lw, linestyle = ls)
#    plt.title("ATTENUATION COEFFICIENT (approximated) \n based on Ainslie & McColm, 1998")
#    plt.legend(title='pH', bbox_to_anchor=(0.2, 0.95), fontsize = 'small')
#    plt.xlabel("Frequency (Hz) \n no pH dependence if S=0!")
#    plt.ylabel("dB km$^{-1}$")
#    plt.tight_layout()
##plt.text(50.3, 60.2, r'$0.106 \frac{f_{1} f^{2}}{f^{2} + f_{1}^{2}} e^{(pH-8)/0.56} + 0.52(1+\frac{T}{43})(\frac{S}{35})\frac{f_{2} f^{2}}{f^{2}+f_{2}^{2}} e^{-z/6}$', horizontalalignment='center', verticalalignment='center', fontsize = 'large', family = 'serif')
##plt.text(50.3, 50.9, r'$+ 0.00049 f^{2} e^{-(T/27 + z/17)}$', horizontalalignment='center', verticalalignment='center', fontsize = 'large', family = 'serif')
#plt.savefig("attenuation_bath_pH.png", dpi = 200)
#plt.show()
#plt.clf()
#
################################################################################
#
## ATTENUATION TEMP- Ainslie & McColm, 1998
#S  = 0. #salinity %. 
#pH = 7.
#T  = [0., 9., 18., 27., 36.] # in #C
#z  = 0.002 # depth in km
## attenuation (dB/km), f in kHz
#
#for i in T:
#    f1 = 0.78*np.sqrt(S/35.)*math.exp(i/26.)
#    f2 = 42.*math.exp(i/17.)
#    lw = 1
#    ls = '--'
#    if i == 18:
#        lw = 2
#        ls = '-'
#    plt.plot(freq, 0.106*(f1*(freq/1000.)**2)/((freq/1000.)**2+f1**2)*math.exp((pH-8.)/0.56) + 0.52*(1+i/43)*(S/35)*(f2*(freq/1000.)**2)/((freq/1000.)**2+f2**2)*math.exp(-z/6) + 0.00049*(freq/1000.)**2*math.exp(-((i/27)+(z/17))), label = str(i), linewidth = lw, linestyle = ls)
#    plt.title("ATTENUATION COEFFICIENT (approximated) \n based on Ainslie & McColm, 1998")
#    plt.legend(title='Temperature', bbox_to_anchor=(0.25, 0.95), fontsize = 'small')
#    plt.xlabel("Frequency (Hz)")
#    plt.ylabel("dB km$^{-1}$")
#    plt.tight_layout()
##plt.text(50.3, 60.2, r'$0.106 \frac{f_{1} f^{2}}{f^{2} + f_{1}^{2}} e^{(pH-8)/0.56} + 0.52(1+\frac{T}{43})(\frac{S}{35})\frac{f_{2} f^{2}}{f^{2}+f_{2}^{2}} e^{-z/6}$', horizontalalignment='center', verticalalignment='center', fontsize = 'large', family = 'serif')
##plt.text(50.3, 50.9, r'$+ 0.00049 f^{2} e^{-(T/27 + z/17)}$', horizontalalignment='center', verticalalignment='center', fontsize = 'large', family = 'serif')
#plt.savefig("attenuation_bath_temp.png", dpi = 200)
#plt.show()
#plt.clf()
#
################################################################################
#
# ATTENUATION DEPTH- Ainslie & McColm, 1998
#S  = 0. #salinity %. 
#pH = 7.
#T  = 18. # in #C
#z  = [0., 0.002, 0.02, 0.2, 2.] # depth in km
#f1 = 0.78*np.sqrt(S/35.)*math.exp(T/26.)
#f2 = 42.*math.exp(T/17.)
## attenuation (dB/km), f in kHz
#
#for i in z:
#    lw = 1
#    ls = '--'
#    if i == 0.002:
#        lw = 2
#        ls = '-'
#    plt.plot(freq, 0.106*(f1*(freq/1000.)**2)/((freq/1000.)**2+f1**2)*math.exp((pH-8.)/0.56) + 0.52*(1+T/43)*(S/35)*(f2*(freq/1000.)**2)/((freq/1000.)**2+f2**2)*math.exp(-i/6) + 0.00049*(freq/1000.)**2*math.exp(-((T/27)+(i/17))), label = str(i), linewidth = lw, linestyle = ls)
#    plt.title("ATTENUATION COEFFICIENT (approximated) \n based on Ainslie & McColm, 1998")
#    plt.legend(title='Depth', bbox_to_anchor=(0.2, 0.95), fontsize = 'small')
#    plt.xlabel("Frequency (Hz)")
#    plt.ylabel("dB km$^{-1}$")
#    plt.tight_layout()
##plt.text(50.3, 60.2, r'$0.106 \frac{f_{1} f^{2}}{f^{2} + f_{1}^{2}} e^{(pH-8)/0.56} + 0.52(1+\frac{T}{43})(\frac{S}{35})\frac{f_{2} f^{2}}{f^{2}+f_{2}^{2}} e^{-z/6}$', horizontalalignment='center', verticalalignment='center', fontsize = 'large', family = 'serif')
##plt.text(50.3, 50.9, r'$+ 0.00049 f^{2} e^{-(T/27 + z/17)}$', horizontalalignment='center', verticalalignment='center', fontsize = 'large', family = 'serif')
#plt.savefig("attenuation_bath_depth.png", dpi = 200)
#plt.show()
#plt.clf()
#
################################################################################
#
#
#########################################
# SITUATIE IN BAD    SITUATIE IN OCEAAN
# ----------------   ----------------
# S  = 0             S  = 38.35
# pH = 7             pH = 7.89
# T  = 18            T  = 13.3
# z  = 0.002         z  = +- 2
##
# https://www.nature.com/articles/srep16770/figures/1
# https://www.nature.com/articles/srep16770
# "Trends of pH decrease in the Mediterranean Sea through high frequency observational data: 
# indication of ocean acidification in the basin", Flecha, S., Pérez, F. F., García-Lafuente3, J., et al., 2015
#########################################
#
################################################################################
#
#plt.clf()
## ATTENUATION - Ainslie & McColm, 1998
#S  = 38.35 #salinity %. 
#pH = 7.89
#T  = 13.3 # in #C
#z  = [0., 1., 2.] # depth in km
#f1 = 0.78*np.sqrt(S/35)*math.exp(T/26)
#f2 = 42.*math.exp(T/17.)
## attenuation (dB/km), f in kHz
#
#for i in z:
#    plt.plot(freq, 0.106*(f1*(freq/1000.)**2)/((freq/1000.)**2+f1**2)*math.exp((pH-8.)/0.56) + 0.52*(1+T/43)*(S/35)*(f2*(freq/1000.)**2)/((freq/1000.)**2+f2**2)*math.exp(-i/6) + 0.00049*(freq/1000.)**2*math.exp(-((T/27)+(i/17))), label = str(int(i))+' km')
#    plt.legend(title='Depth (z)', bbox_to_anchor=(0.97, 0.35), fontsize = 'small')
#    plt.title("ATTENUATION COEFFICIENT (approximated) \n based on Ainslie & McColm, 1998")
#    plt.xlabel("Frequency (Hz)")
#    plt.ylabel("dB km$^{-1}$")
#plt.tight_layout()
#plt.text(0.4, 1.84,'Sea properties:\n S = %s '% S + '\n pH = %s'% pH + '\n T = %s '% T, transform=ax.transAxes, horizontalalignment='left', verticalalignment='center', fontsize = 'small', family = 'serif', bbox = dict(facecolor='white', edgecolor='black', pad=10.0))
#
##plt.text(10.3, 2.2, r'$0.106 \frac{f_{1} f^{2}}{f^{2} + f_{1}^{2}} e^{(pH-8)/0.56} + 0.52(1+\frac{T}{43})(\frac{S}{35})\frac{f_{2} f^{2}}{f^{2}+f_{2}^{2}} e^{-z/6}$', horizontalalignment='center', verticalalignment='center', fontsize = 'large', family = 'serif')
##plt.text(3000, 30, r'$+ 0.00049 f^{2} e^{-(T/27 + z/17)}$', horizontalalignment='center', verticalalignment='center', fontsize = 'large', family = 'serif')
#plt.savefig("attenuation_sea.png", dpi = 200)
#plt.show()
#plt.clf()
#

###########################################
###########################################
#####                                 #####
##### - Kurahashi, N., et al., 2008 - #####
#####                                 #####
###########################################
###########################################

################################################################################

# - AMBIENT NOISE AT SS0
plt.plot(freq, np.log10((freq/1000.)**(-5./3.)) + 94.5, label='10log')
#plt.legend(title='seastate0', bbox_to_anchor=(1.04, 1.05), fontsize = 'small')
plt.title("AMBIENT NOISE LEVEL FOR SEASTATE 0 \n based on Kurahashi et al., 2008")
plt.xlabel("Frequency (Hz)")
plt.ylabel("dB re 1 $\mu Pa^{2}/Hz$")
plt.xlim([1000,15000])
plt.ylim([92, 95])
plt.tight_layout()
#plt.text(1.3, 1.3,'NL = log10 (f$^{-5/3}$) +94.5', horizontalalignment='center', verticalalignment='center', transform=ax.transAxes, fontsize = 'medium', family = 'serif')
plt.savefig("seastate0_kurahashi_2008.png", dpi=200)
plt.show()
plt.clf()
#
# - AMBIENT NOISE AT SS0 sqrt taken
#
#plt.plot(freq, np.sqrt(np.log10((freq/1000.)**(-5./3.)) + 94.5), label='10log')
#plt.legend(title='seastate0', bbox_to_anchor=(1.04, 1.05), fontsize = 'small')
#plt.title("AMBIENT NOISE LEVEL SEASTATE 0 \n based on Kurahashi et al., 2008")
#plt.xlabel("Frequency (Hz)")
#plt.ylabel("$\sqrt{dB re 1 \mu Pa^{2}/Hz}$")
##plt.xlim([1000,25000])
##plt.ylim([30, 70])
#plt.tight_layout()
#plt.text(1.3, 1.3,'NL = $\sqrt{log10 (f^{-5/3}) +94.5}$ \n ik denk dat je niet zomaar de wortel mag nemen.', horizontalalignment='center', verticalalignment='center', transform=ax.transAxes, fontsize = 'medium', family = 'serif')
#plt.savefig("seastate0_sqrt_kurahashi_2008.png", dpi=200)
#plt.show()
#plt.clf()
#
# - SEASTATE AMBIENT NOISE
#
ss       = [0., 1., 2., 3., 4., 5., 6.]
for i in ss:
    plt.plot(freq, np.log10((freq/1000.)**(-5./3.)) + 94.5 + 30 * np.log10(i + 1), label=str(int(i)))
    plt.legend(title='seastate', bbox_to_anchor=(1.04, 1.05), fontsize = 'small')
    plt.title("AMBIENT NOISE LEVEL PER SEASTATE \n based on Kurahashi et al., 2008")
    plt.xlabel("Frequency (Hz)")
    plt.ylabel("dB re 1 $\mu Pa^{2}/Hz$")
    plt.xlim([1000,15000])
#    plt.ylim([30, 70])
    plt.tight_layout()
#plt.text(1.3, 1.8,'NL = 56 + 19 log (ss) - 17 log (f)', horizontalalignment='center', verticalalignment='center', fontsize = 'medium', family = 'serif')
plt.savefig("seastates_Kurahashi_1984.png", dpi=200)
plt.show()
plt.clf()
