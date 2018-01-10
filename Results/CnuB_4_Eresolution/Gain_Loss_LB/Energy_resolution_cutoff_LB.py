# -*- coding: utf-8 -*-
"""
Created on Thu Nov 16 15:59:51 2017

@author: mullerrs
"""
import os
import sys
from math import *
import numpy as np
import matplotlib.pyplot as plt
import ROOT
ROOT.gROOT.SetBatch(True)   # so it doesn't try to actively popup graphics on the screen <= needed for stoomboot
from scipy.interpolate import interp1d
import scipy.integrate as integrate
from bisect import bisect_left

filename = os.path.splitext(os.path.basename(sys.argv[1]))[0]
basename = filename.split("H0H1")[0]

# =================
# number of events
# =================

def events(Flux):
    '''Function to convert flux to amount of events'''
    Detector_eff = 0.01 # [fraction] for now assume 1% detector efficiency       <== Gauss
    t = 10              # [years]
    A = pi*500**2       # [m^2]      pi * (500m)^2 = opp 1 Antares blok
    Theta = 2*pi        # [sr]       Full sphere is 4pi sr
                        #            --> here half of the angular resolution
    return (60*60*24*365.25 * t) * A * Theta * Detector_eff * Flux


# =================
# Flux data
# =================

def spectrum_func(e, H, E_H, Flux_H0, Flux_H1 ):    # mean expected events
    '''Function that returns a function for the null hypothesis H = H0 
    (constant spectrum), and H = H1 hypothesis (dip in spectrum) '''

    # interpolate the separate datapoints
    f_H0 = interp1d(E_H, Flux_H0)
    f_H1 = interp1d(E_H, Flux_H1)

    if H == 'H0':
        return f_H0(e)
    elif H == 'H1':
        return f_H1(e)


def fill_hist( func, H, E_H, Flux_H0, Flux_H1):
    '''Function to plot a histogram according to a function func, given the 
    hypothesis H'''
    bins = 100
    h1 = ROOT.TH1D( 'h1','Flux', bins, min(E_H), max(E_H))

    Binwidth = (max(E_H)-min(E_H))/bins

    for be in range( 1, h1.GetNbinsX() + 1 ):
        e = h1.GetXaxis().GetBinCenter( be )
        b = h1.GetBin( be ) # global bin number
        z = func(e, H, E_H, Flux_H0, Flux_H1 )
        h1.SetBinContent( b, z )

    h1.SetXTitle("Energy (eV)")
    h1.SetYTitle("Flux (Arbitrary unit)")

    return h1


def draw_H0_H1(what, E_H, Flux_H0, Flux_H1):
    '''Function that draws H0 and H1 from imported data in a histogram'''
    c1 = ROOT.TCanvas()
    h1 = fill_hist( spectrum_func , 'H0', E_H, Flux_H0, Flux_H1)
    h1.SetLineStyle(3)
    h1.Draw('LSAME')
    h1.SetTitle("Flux function")
    h2 = fill_hist( spectrum_func , 'H1', E_H, Flux_H0, Flux_H1)
    h2.Draw('LSAME')
    h2.SetTitle("Flux function")
    c1.Update()

    h1.Write("Original_Flux_H0"+what)
    h2.Write("Original_Flux_H1"+what)
    c1.Write("Original_Flux_canvas"+what)

    h1.Delete()
    h2.Delete()

# =================
# Inverse Transform
# =================

def binary_search(a, x):
    '''Function that returns the position in list a closest to value x
    if x is exactly halfway two values, it takes the upper index,
    if x is below the lowest value, index 0 is returned
    if x is above the upper value, the last index is returned'''
    if x < a[0]:
        return 0
    elif x > a[-1]:
        return -1
    else:
        pos = bisect_left(a, x)  # find insertion position coming from left
        if float(a[pos]) == float(x):
            return pos
        elif a[pos] - x > x - a[pos-1]:
            return pos - 1
        else:
            return pos


def integr_spectrum_func(E_H, Flux_H):
    '''Function thet returns stepwise integrated Flux'''
    Int_Flux_H = [Flux_H[0]]
    for i in range(len(E_H)-1):
        i += 1
        Int_Flux_H.append(Int_Flux_H[-1]+Flux_H[i])
    return Int_Flux_H


def inverse_transform(E_H, Int_Flux_H0, Int_Flux_H1, H):
    if H == 'H1': 
        Int_Flux = Int_Flux_H1
    elif H == 'H0':
        Int_Flux = Int_Flux_H0

    # bounds of the fluxvalues to create random number in:
    Int_Flux_max = max(Int_Flux)
    
    rand_f = ROOT.gRandom.Rndm() * Int_Flux_max                                   # uniformly distributed number between 0 and Flux_max
    rand_e = E_H[ binary_search(Int_Flux, rand_f) ]                               # corresponding closest defined Energy

    return rand_e

# =================
# Likelihood
# =================

def log_likelihood(hist, H, N_events, E_H, Flux_H0, Flux_H1):
    '''Function that returns the log likelihood for a given hypothesis 'H' with
    respect to the given histogram 'hist' with your data. 
    NB: make sure Flux_H0 and Flux_H1, interpolated are normalized to 1'''
    loglik = 0.

    # -- loop over bins
    for i_bin in range( 1, hist.GetNbinsX() + 1 ):
        e = hist.GetBinCenter( i_bin )                                          # energy (centre of bin) 
        mu_bin = spectrum_func(e, H, E_H, Flux_H0, Flux_H1) * hist.GetBinWidth(1) * N_events     # theoretical amount of counts in bin
        Nevt_bin = hist.GetBinContent( i_bin )                                  # the amount of counts in bin found in the data of histogram h

        if ROOT.TMath.Poisson( Nevt_bin, mu_bin ) > 0:                          # <= check if&when this will happen
            loglik += ROOT.TMath.Log( ROOT.TMath.Poisson( Nevt_bin, mu_bin ) )
        else:
            print "Negative probability!!"
    return loglik


def LLR(hist, N_events, E_H, Flux_H0, Flux_H1):
    '''Function that caluclates the likelihood ratio'''
    L_H0 = log_likelihood(hist, 'H0', N_events, E_H, Flux_H0, Flux_H1)
    L_H1 = log_likelihood(hist, 'H1', N_events, E_H, Flux_H0, Flux_H1)
    LLR = L_H1 - L_H0
    return LLR


def plot_LLR_value_in_hist(N_events, bins, Eresolution, H, hist, E_H, Flux_H0, Flux_H1, Int_Flux_H0, Int_Flux_H1):
    '''Function that fills the given histogram 'hist', with the LLR for 
    pseudo experiment data based on hypothesis 'H' '''
    h4 = pseudo_exp(N_events, bins, Eresolution, H, E_H, Int_Flux_H0, Int_Flux_H1)
    LLratio_Hdata = LLR(h4, N_events, E_H, Flux_H0, Flux_H1)
    hist.Fill(LLratio_Hdata)
    h4.Delete()
    return LLratio_Hdata


def pseudo_exp(N_events, bins, Eresolution, H, E_H, Int_Flux_H0, Int_Flux_H1):
    '''Functin that creates pseudo experiments of N_events detections,
    based on the acceptence - rejection method'''
#    c3 = ROOT.TCanvas()
    h3 = ROOT.TH1D( 'h3','Pseudo Experiment, N=%s'%(N_events), bins, min(E_H), max(E_H))
    for i in range(N_events):
#        if i % (N_events/10) == 0:
#            print "%s/%s" %(i, N_events)
        E = inverse_transform(E_H, Int_Flux_H0, Int_Flux_H1, H)
        E = E * ROOT.gRandom.Gaus(1., Eresolution/100.) #Gauss(mean, sigma)          # <= wat doen met E<0?!
        h3.Fill(E)
#    h3.Draw()
    h3.SetXTitle("Energy (eV)")
    h3.SetYTitle("Counts per bin")
    
#    c3.SaveAs("inverse_transform_pseudo_event_Eresol%s_N%s_%s.png"%(Eresolution, N_events, H))

    h3.Write('Pseudo_exp_InvTr, %s, N=%s, Eres=%s'%(H, N_events, Eresolution))
    return h3


def determine_p_value(h_h0, teststatistic_data):
    '''Funtion to determine the p-value of the data with given test statistic: 
    determine P(teststatistic > teststatistic_data | H0) by integration'''
    # integrate from tets_statistic to last bin
    axis      = h_h0.GetXaxis()
    bin_min   = axis.FindBin(teststatistic_data)
    bin_max   = h_h0.GetNbinsX() + 1
    integral  = h_h0.Integral(bin_min, bin_max)
    # subtract part of the first bin that was too much
    integral -= h_h0.GetBinContent(bin_min) * ( teststatistic_data - axis.GetBinLowEdge(bin_min) ) / axis.GetBinWidth(bin_min)
    # normalize
    p_value = integral/h_h0.GetEntries()
    return p_value


def determine_CL_value(h_h1, teststatistic_data):
    '''Funtion to determine the CL-value of the data with given test statistic: 
    determine P(teststatistic < teststatistic_data | H1) by integration'''
    # integrate from first bin to tets_statistic
    axis      = h_h1.GetXaxis()
    bin_min   = 1
    bin_max   = axis.FindBin(teststatistic_data)
    integral  = h_h1.Integral(bin_min, bin_max)
    # subtract part of the last bin that was too much
    integral -= h_h1.GetBinContent(bin_max) * ( (axis.GetBinLowEdge(bin_max) + h_h1.GetBinWidth(bin_max)) - teststatistic_data ) / axis.GetBinWidth(bin_min)
    # normalize
    CL_value = integral/h_h1.GetEntries()
    return CL_value


def plot_Hypothesis_test(h_h0, h_h1, savingname):
    cc = ROOT.TCanvas()

    #   Make up the plot
    h_h0.Draw()
    h_h0.SetLineColor(1)
    h_h1.Draw('SAME')
    h_h1.SetLineColor(3)
    #cc.BuildLegend()
    h_h1.SetTitle('Hypothesis Testing')
    h_h0.SetTitle('Hypothesis Testing')
    h_h0.SetXTitle("Test statistic log(lambda)")
    h_h0.SetYTitle("Probability density (counts)")

    h_h0.Write('Htest_%s_N=%s_Eres=%s'%('H0', N_events, Eresolution))
    h_h1.Write('Htest_%s_N=%s_Eres=%s'%('H1', N_events, Eresolution))

    cc.SaveAs(savingname)


def plot_line(canvas, hist, teststatistic, savingname):
    '''Function that plots the line of the test statistic of the pseudo experiment'''
    line = ROOT.TLine(teststatistic, 0, teststatistic, 3/5.*hist.GetMaximum())
    line.SetLineColor(2)                                                            # for shaded area: https://root.cern.ch/root/roottalk/roottalk98/0846.html
    line.SetLineWidth(3)
    line.Draw('SAME')
    canvas.Update()
    canvas.SaveAs(savingname)


def median(lst):
    '''Function to determine the median value of a given list 'lst' '''
    quotient, remainder = divmod(len(lst), 2)
    if remainder:
        return sorted(lst)[quotient]
    return sum(sorted(lst)[quotient - 1:quotient + 1]) / 2.


# =================
# main & usage
# =================

def usage():
    print "Usage:  python  %s  <Fluxdata-file> \n" % os.path.basename(sys.argv[0])

def main():

    f = ROOT.TFile(basename+"_histos.root", "recreate")                         # create root file to save everything in

    ROOT.gStyle.SetOptStat(0)                                                   # do not show statistics box
    ROOT.gStyle.SetTitleOffset(1.3, "y")                                        # spacing y-label

    args = sys.argv[1:]
    if "-h" in args or "--help" in args or len(args) < 1:
        usage()
        sys.exit(2)

    # -- Import flux data
    E_H, Flux_H0, Flux_H1 = np.loadtxt(sys.argv[1], delimiter = '\t', unpack=True)
    draw_H0_H1("Flux", E_H, Flux_H0, Flux_H1)

    # -- Normalize flux data (and plot)
    Binsize = E_H[1]-E_H[0]
    Flux_H0_norm = Flux_H0 / (sum(Flux_H0) * Binsize)                           # Also divide by binsize to make sure that the interpolated function is normalized to 1
    Flux_H1_norm = Flux_H1 / (sum(Flux_H1) * Binsize)
    draw_H0_H1("Norm", E_H, Flux_H0_norm, Flux_H1_norm)

    # -- Check that the continuous function made with spectrum_func given Flux_H_norm data equals 1
#    f0_to_integrate = lambda e: spectrum_func(e, 'H0', E_H, Flux_H0_norm, Flux_H1_norm)
#    f1_to_integrate = lambda e: spectrum_func(e, 'H1', E_H, Flux_H0_norm, Flux_H1_norm)
#    integrant0, err1 = integrate.quad(f0_to_integrate, min(E_H), max(E_H))
#    integrant1, err1 = integrate.quad(f1_to_integrate, min(E_H), max(E_H))
#    print 'integrant0, integrant1', integrant0, integrant1

    # -- Integrate flux data for inverse transform method
    Int_Flux_H0 = integr_spectrum_func(E_H, Flux_H0_norm)
    Int_Flux_H1 = integr_spectrum_func(E_H, Flux_H1_norm)
    draw_H0_H1("Int", E_H, Int_Flux_H0, Int_Flux_H1)

    N_events = 100000
    # -- Create Pseudo Measurement for N_events according to H0, and H1, (plot, and save figure)
    print 'test 1'
    pseudo_exp(N_events, 100, 0, 'H0', E_H, Int_Flux_H0, Int_Flux_H1)
    pseudo_exp(N_events, 100, 0, 'H1', E_H, Int_Flux_H0, Int_Flux_H1)
    print 'test 2'
    pseudo_exp(N_events, 100, 30, 'H0', E_H, Int_Flux_H0, Int_Flux_H1)
    pseudo_exp(N_events, 100, 30, 'H1', E_H, Int_Flux_H0, Int_Flux_H1)
    print 'test 3'
    pseudo_exp(N_events, 100, 80, 'H0', E_H, Int_Flux_H0, Int_Flux_H1)
    pseudo_exp(N_events, 100, 80, 'H1', E_H, Int_Flux_H0, Int_Flux_H1)


    # =================
    # -- Perform Eresolution test
    # =================

    Elb = 5e20
    Eub = 5e22
#    N_events = int(N_events*(max(E_H)-min(E_H))/(Eub-Elb))
    bins = int(10000*(max(E_H)-min(E_H))/(Eub-Elb))
    EResolution = [0., 30., 80.]
    H = 'H1'
    
    for Eresolution in EResolution:
        print Eresolution
        canvas = ROOT.TCanvas()
    
        ''' Based on function "pseudo experiment" '''
        HistOud = ROOT.TH1D( 'HistOud','Pseudo Experiment, N=%s'%(N_events), bins, min(E_H), max(E_H))
        HistNieuw = ROOT.TH1D( 'HistNieuw','Pseudo Experiment, N=%s'%(N_events), bins, min(E_H), max(E_H))
        LinksBij_oud = ROOT.TH1D( 'LinksBij_oud','Pseudo Experiment, N=%s'%(N_events), bins, min(E_H), max(E_H))
        LinksBij_nieuw = ROOT.TH1D( 'LinksBij_nieuw','Pseudo Experiment, N=%s'%(N_events), bins, min(E_H), max(E_H))
        LinksWeg_oud = ROOT.TH1D( 'LinksWeg_oud','Pseudo Experiment, N=%s'%(N_events), bins, min(E_H), max(E_H))
        LinksWeg_nieuw = ROOT.TH1D( 'LinksWeg_nieuw','Pseudo Experiment, N=%s'%(N_events), bins, min(E_H), max(E_H))
        RechtsBij_oud = ROOT.TH1D( 'RechtsBij_oud','Pseudo Experiment, N=%s'%(N_events), bins, min(E_H), max(E_H))
        RechtsBij_nieuw = ROOT.TH1D( 'RechtsBij_nieuw','Pseudo Experiment, N=%s'%(N_events), bins, min(E_H), max(E_H))
        RechtsWeg_oud = ROOT.TH1D( 'RechtsWeg_oud','Pseudo Experiment, N=%s'%(N_events), bins, min(E_H), max(E_H))
        RechtsWeg_nieuw = ROOT.TH1D( 'RechtsWeg_nieuw','Pseudo Experiment, N=%s'%(N_events), bins, min(E_H), max(E_H))
    
        Hist_Helemaal_oud = ROOT.TH1D( 'Hist_Helemaal_oud','Pseudo Experiment, N=%s'%(N_events), bins, min(E_H), max(E_H))
        Hist_Links_oud = ROOT.TH1D( 'Hist_Links_oud','Pseudo Experiment, N=%s'%(N_events), bins, min(E_H), max(E_H))
        Hist_Rechts_oud = ROOT.TH1D( 'Hist_Rechts_oud','Pseudo Experiment, N=%s'%(N_events), bins, min(E_H), max(E_H))

        for i in range(N_events):
            if i % (N_events/100) == 0:
                print "%s/%s" %(i, N_events)
    
            Eoud = inverse_transform(E_H, Int_Flux_H0, Int_Flux_H1, H)
            Enieuw = Eoud * ROOT.gRandom.Gaus(1., Eresolution/100.) #Gauss(mean, sigma)          # <= wat doen met E<0?!

            Hist_Helemaal_oud.Fill(Eoud)
            if Eoud < Elb:
                Hist_Links_oud.Fill(Eoud)
            if Eoud > Eub:
                Hist_Rechts_oud.Fill(Eoud)
    
            # Wat was de oorspronkelijke histogram, en wordt de nieuwe? MID ONLY
            if Eoud >= Elb and Eoud <= Eub:
                HistOud.Fill(Eoud)
                if Enieuw >= Elb and Enieuw <= Eub:
                    HistNieuw.Fill(Enieuw)
    
            # Welke komen er van links bij?
            if Eoud < Elb and Enieuw >= Elb and Enieuw <= Eub:
                LinksBij_oud.Fill(Eoud)         #licht blauw
                LinksBij_nieuw.Fill(Enieuw)     #donker blauw
            # Welke verlies je links?
            elif Eoud >= Elb and Eoud <= Eub and Enieuw < Elb:
                LinksWeg_oud.Fill(Eoud)         # licht rood
                LinksWeg_nieuw.Fill(Enieuw)     # donker rood
            # Welke komen er rechts bij?
            elif Eoud > Eub and Enieuw <= Eub and Enieuw >= Elb:
                RechtsBij_oud.Fill(Eoud)        # licht geel
                RechtsBij_nieuw.Fill(Enieuw)    # donker geel
            # Welke verlies je rechts?
            elif Eoud <= Eub and Eoud >= Elb and Enieuw > Eub:
                RechtsWeg_oud.Fill(Eoud)        # licht groen
                RechtsWeg_nieuw.Fill(Enieuw)     # donker groen
    
        LinksBij_oud.SetLineColor(600)
        LinksBij_nieuw.SetFillColor(600)
        LinksBij_nieuw.SetLineColor(600)
        LinksWeg_oud.SetLineColor(632)
        LinksWeg_nieuw.SetFillColor(632)
        LinksWeg_nieuw.SetLineColor(632)
        RechtsBij_oud.SetLineColor(400)
        RechtsBij_nieuw.SetFillColor(400)
        RechtsBij_nieuw.SetLineColor(400)
        RechtsWeg_oud.SetLineColor(416)
        RechtsWeg_nieuw.SetFillColor(416)
        RechtsWeg_nieuw.SetLineColor(416)
        HistOud.SetLineColor(1) #zwart
        HistNieuw.SetLineColor(920) #grijs
    
        HistOud.SetXTitle("Energy (eV)")
        HistOud.SetYTitle("Counts per bin")
    
        HistOud.Write('HistOud, %s, N=%s, Eres=%s'%(H, N_events, Eresolution))
        HistNieuw.Write('HistNieuw, %s, N=%s, Eres=%s'%(H, N_events, Eresolution))
        LinksBij_oud.Write('LinksBij_oud, %s, N=%s, Eres=%s'%(H, N_events, Eresolution))
        LinksBij_nieuw.Write('LinksBij_nieuw, %s, N=%s, Eres=%s'%(H, N_events, Eresolution))
        LinksWeg_oud.Write('LinksWeg_oud, %s, N=%s, Eres=%s'%(H, N_events, Eresolution))
        LinksWeg_nieuw.Write('LinksWeg_nieuw, %s, N=%s, Eres=%s'%(H, N_events, Eresolution))
        RechtsBij_oud.Write('RechtsBij_oud, %s, N=%s, Eres=%s'%(H, N_events, Eresolution))
        RechtsBij_nieuw.Write('RechtsBij_nieuw, %s, N=%s, Eres=%s'%(H, N_events, Eresolution))
        RechtsWeg_oud.Write('RechtsWeg_oud, %s, N=%s, Eres=%s'%(H, N_events, Eresolution))
        RechtsWeg_nieuw.Write('RechtsWeg_nieuw, %s, N=%s, Eres=%s'%(H, N_events, Eresolution))
    
        HistOud.Draw()
        HistNieuw.Draw('SAME')
        LinksBij_oud.Draw('SAME')
        LinksBij_nieuw.Draw('SAME')
        LinksWeg_oud.Draw('SAME')
        LinksWeg_nieuw.Draw('SAME')
        RechtsBij_oud.Draw('SAME')
        RechtsBij_nieuw.Draw('SAME')
        RechtsWeg_oud.Draw('SAME')
        RechtsWeg_nieuw.Draw('SAME')
        
        linelb = ROOT.TLine(Elb, 0, Elb, HistOud.GetMaximum())
        linelb.SetLineColor(920)
        linelb.SetLineWidth(3)
        linelb.Draw('SAME')
    
        lineub = ROOT.TLine(Eub, 0, Eub, HistOud.GetMaximum())
        lineub.SetLineColor(920)
        lineub.SetLineWidth(3)
        lineub.Draw('SAME')

        canvas.Update()
        canvas.Write("Analyze_Eresol%s_smeering_N%s_%s"%(Eresolution, N_events, H))

        canvas2 = ROOT.TCanvas()
        Hist_Helemaal_oud.Draw()
        Hist_Helemaal_oud.SetLineColor(1)
        LinksBij_nieuw.Draw('SAME')
        LinksBij_oud.Draw('SAME')
        linelb.Draw('SAME')
        canvas.Update()
        canvas2.Write("Analyze_Eresol%s_smeering_N%s_%s_LINKS_BIJ"%(Eresolution, N_events, H))

        canvas3 = ROOT.TCanvas()
        Hist_Helemaal_oud.Draw()
        Hist_Helemaal_oud.SetLineColor(1)
        LinksWeg_nieuw.Draw('SAME')
        LinksWeg_oud.Draw('SAME')
        linelb.Draw('SAME')
        canvas.Update()
        canvas3.Write("Analyze_Eresol%s_smeering_N%s_%s_LINKS_WEG"%(Eresolution, N_events, H))

        print "Eresolution = %s" %(Eresolution)

        print "\nLowerbound = %i" %Elb
        print "Number of events below lower bound = %i, compared to the total of %i thrown events." %(Hist_Links_oud.GetEntries(), Hist_Helemaal_oud.GetEntries())
        print "Number of events from above the lower bound, fallen below this bound \n= ", LinksWeg_oud.GetEntries(), " = ",  LinksWeg_nieuw.GetEntries()
        print "Number of events from below the lower bound, fallen above this bound \n= ", LinksBij_oud.GetEntries(), " = ", LinksBij_nieuw.GetEntries()
    
        print "\nUpperbound = %i" %Eub
        print "Number of events above upper bound = %i, compared to the total of %i thrown events." %(Hist_Rechts_oud.GetEntries(), Hist_Helemaal_oud.GetEntries())
        print "Number of events from below the upeer bound, fallen above this bound \n= ", RechtsWeg_oud.GetEntries(), " = ", RechtsWeg_nieuw.GetEntries()
        print "Number of events from above the upper bound, fallen below this bound \n= ", RechtsBij_oud.GetEntries(), " = ", RechtsBij_nieuw.GetEntries()

    f.Close()

if __name__ == "__main__":
    sys.exit(main())

