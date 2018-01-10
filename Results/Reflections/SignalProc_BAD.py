"""

This module contains the Signal processing class for estimating
non-parametric models

"""
import math
from numpy.linalg import norm
from numpy.fft import fft
from numpy import cos
from numpy import pi
import numpy as np
from scipy import signal

import matplotlib.pyplot as plt
class SignalProc(object):

    """

    A class for performing non-parametric model fitting

    .. note::   Ported to Python v3.x.x
                and therfore not more downwards compatible with Python v 2.x.x

    Parameters
    ----------
    Nchnl : int
        The number of signals.
    length : int
        The size of the FFT segment.
    Fs : float
        The sampling frequency of the acquired time series
    
    Methods
    -------
    MakeSig(self, F0):
        Generate a sin(2*pi*F0/Fs*time)
    .. math:: X(e^{j\omega } ) = x(n)e^{ - j\omega n}
    Spectrum()
        Caluculate the Power Density Spectrum

    Returns
    -------
    An Instance of the class

    Examples
    --------
    >>> ss = SignalProc(6,10240,2048)
    >>> ss.MakeSig(128)
    >>> ss.WindowType ='hanning'
    >>> ss.NFFT       = 1024
    >>> ss.PSDOn      = 1
    >>> ss.Spectrum(12)
    >>> ss.plotSpectrum()

    Copyright
    ---------
    E.J.J. Doppenberg
    IPOS.vision
    12-05-2016

    """

    Counter = 0;
    def __init__(self, *args):
        self.basename   = 'basename'
        self.refr       = 1 
        self.resp       = 0
        nargin          = len(args)
        self.F0         = 63.
        self.norm       = 1.
        self.NFFT       = 1024
        self.Segment    = 100
        self.PSDOn      = 0
        self.OnePlot    = 0
        self.Cummulative= 0
        self.Fs         = 1024.
        self.length     = 102400.
        self.Nchnl      = 4
        self.FcumBegin  = 0
        self.FcumEnd    = self.Fs/2
        if nargin >= 1:
            self.Nchnl  = args[0]
        if nargin >= 2:
            self.length = args[1]
        if nargin >= 3:
            self.Fs      = args[2]
        self.WindowType  = 'hanning'
        self.window      = np.hanning(self.NFFT)

        self.CalData     = np.zeros((self.length, self.Nchnl))
        self.Data        = np.zeros((self.length, self.Nchnl))
        self.freq        = self.MakeFreqSpan()
        self.time        = self.MakeTimeSpan()
        self.PSD         = np.ones((self.NFFT/2, self.Nchnl), dtype=float)
        self.Sx          = np.zeros((self.NFFT, self.Nchnl), dtype=complex)
        self.Sy          = np.zeros((self.NFFT, self.Nchnl), dtype=complex)
        self.Sxy         = np.zeros((self.NFFT, self.Nchnl), dtype=complex)
        self.Txy         = np.zeros((self.NFFT, self.Nchnl), dtype=complex)
        self.Sxx         = np.zeros((self.NFFT, self.Nchnl), dtype=float)
        self.Syy         = np.zeros((self.NFFT, self.Nchnl), dtype=float)
        self.Cxy         = np.zeros((self.NFFT, self.Nchnl), dtype=float)
        self.gg          = np.ones((self.Nchnl), dtype=float)
        self.Gain        = np.ones((self.Nchnl), dtype=float)
        self.Sens        = np.ones((self.Nchnl), dtype=float)
        self.Unit        = ['V']
        self.Name        = ['P']
        self.GraphLabel  = ['Exp']
        self.Title       = ['Graph']
        for nchl in range(1,self.Nchnl):
            self.Unit       += ['V']
            self.Name       += ['P']
            self.GraphLabel += ['Exp']
            self.Title      += ['Graph']

        SignalProc.Counter +=1
 #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    def __call__(self):
        self.NFFT       = 4096
        self.Segment    = 100
        self.Fs         = 4096
        self.Nchnl      = 4
        self.FcumBegin  = 0
        self.FcumEnd    = self.Fs/2
        self.WindowType = 'hanning'
        self.window     = self.MakeWindow()
        self.length     = 102400
        self.PSDOn      = 0
        self.OnePlot    = 0
        self.Cummulative= 0
        self.CalData    = np.zeros(self.length, self.Nchnl)
        self.Data       = np.zeros(self.length, self.Nchnl)
        self.freq       = self.MakeFreqSpan()
        self.time       = self.MakeTimeSpan()
        self.PSD        = np.ones((self.NFFT/2, self.Nchnl), dtype=float)
        self.Sx         = np.zeros((self.NFFT, self.Nchnl), dtype=complex)
        self.Sy         = np.zeros((self.NFFT, self.Nchnl), dtype=complex)
        self.Sxy        = np.zeros((self.NFFT, self.Nchnl), dtype=complex)
        self.Txy        = np.zeros((self.NFFT, self.Nchnl), dtype=complex)
        self.Sxx        = np.zeros((self.NFFT, self.Nchnl), dtype=float)
        self.Syy        = np.zeros((self.NFFT, self.Nchnl), dtype=float)
        self.Cxy        = np.zeros((self.NFFT, self.Nchnl), dtype=float)
        self.Gain        = np.ones((self.Nchnl), dtype=float)
        self.Sens        = np.ones((self.Nchnl), dtype=float)
        self.Unit        = ['V']
        self.Name        = ['P']
        self.GraphLabel  = ['Exp']
        self.Title       = ['Graph']
        for nchl in range(1, self.Nchnl):
            self.Unit       += ['V']
            self.Name       += ['P']
            self.GraphLabel += ['Exp']
            self.Title      += ['Graph']

        return  # dit doet toch niks?
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    def __repr__(self):
        Rspace = 30;
        str3  ="Fs:".rjust(Rspace)+" %d\n".ljust(0) %self.Fs
        str3 += "NFFT:".rjust(Rspace)+" %d\n".ljust(0) %self.NFFT
        str3 += "Segment:".rjust(Rspace)+" %d\n".ljust(0) %self.Segment
        str3 += "PSDOn:".rjust(Rspace)+" %d\n".ljust(0) %self.PSDOn
        str3 += "OnePlot:".rjust(Rspace)+" %d\n".ljust(0) %self.OnePlot
        str3 += "Cummulative:".rjust(Rspace)+" %d\n".ljust(0) %self.Cummulative
        str3 += "length:".rjust(Rspace)+" %d\n".ljust(0) %self.length
        str3 += "Nchnl:".rjust(Rspace)+" %d\n".ljust(0) %self.Nchnl
        str3 += "FCumBegin:".rjust(Rspace)+" %d\n".ljust(0) %self.FcumBegin
        str3 += "FCumEnd:".rjust(Rspace)+" %d\n".ljust(0) %self.FcumEnd
        str3 += "refr:".rjust(Rspace)+" %d\n".ljust(0) %self.refr
        str3 += "resp:".rjust(Rspace)+" %d\n".ljust(0) %self.resp
        str3 += "WindowType:".rjust(Rspace)+" '%s'\n".ljust(0) %self.WindowType
        str3 += "window:".rjust(Rspace)+" [%d]\n".ljust(0) % np.size(self.window)
        str3 += "CalData:".rjust(Rspace)+" [%s]\n".ljust(0) % str(np.shape(self.CalData))
        str3 += "Data:".rjust(Rspace)+" [%s]\n".ljust(0) % str(np.shape(self.Data))
        str3 += "time:".rjust(Rspace)+" [%d]\n".ljust(0) % np.size(self.time)
        str3 += "freq:".rjust(Rspace)+" [%d]\n".ljust(0) % np.size(self.freq)
        str3 += "PSD:".rjust(Rspace)+" [%s]\n".ljust(0) %str(np.shape(self.PSD))
        str3 += "Sx:".rjust(Rspace)+" [%s]\n".ljust(0) % str(np.shape(self.Sx))
        str3 += "Sy:".rjust(Rspace)+" [%s]\n".ljust(0) % str(np.shape(self.Sy))
        str3 += "Sxx:".rjust(Rspace)+" [%s]\n".ljust(0) % str(np.shape(self.Sx))
        str3 += "Syy:".rjust(Rspace)+" [%s]\n".ljust(0) % str(np.shape(self.Sy))
        str3 += "Sxy:".rjust(Rspace)+" [%s]\n".ljust(0) % str(np.shape(self.Sx))
        str3 += "Txy:".rjust(Rspace)+" [%s]\n".ljust(0) % str(np.shape(self.Sy))
        str3 += "Cxy:".rjust(Rspace)+" [%s]\n".ljust(0) % str(np.shape(self.Sx))
        str3 += "GraphLabel:".rjust(Rspace)+" [%s]\n".ljust(0) % self.GraphLabel
        str3 += "Unit:".rjust(Rspace)+" [%s]\n".ljust(0) % self.Unit
        str3 += "Name:".rjust(Rspace)+" [%s]\n".ljust(0) % self.Name
        str3 += "Title:".rjust(Rspace)+" [%s]\n".ljust(0) % self.Title
        str3 += "Gain:".rjust(Rspace)+"[%s]\n".ljust(0) % str(self.Gain)
        str3 += "Sens:".rjust(Rspace)+"[%s]\n".ljust(0) % str(self.Sens)

        return str3;
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    def Calibrate(self):
#      self.CalData = np.copy(self.Data)
      self.gg           = 1./(self.Gain*self.Sens)
      self.CalData = np.multiply(self.gg, self.Data)
      return self.CalData
 #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    def plot(self, rChnl):
      self.Calibrate()
      self.length = np.size(self.CalData,0)
      self.MakeTimeSpan()
      if self.OnePlot == 0:
          nPlot    = len(rChnl)
          nSubplot = 1
          for nChnl in rChnl:
              plt.subplot(nPlot, 1, nSubplot)
              nSubplot += 1

#              ##################
#              #### alle data ####
#              ##################
#              plt.plot(self.time, self.CalData[:,nChnl])
#
#              ##################
#              #### 1 sec ####
#              ##################
#              plt.plot(self.time[:self.Fs], self.CalData[:self.Fs,nChnl]) # <-- 1 window plotten
#              plt.xlabel('#samples, \n NB: 1 sample = 5.12e-6 sec') # <<- x-axis for plotting 1 window
#              # NB: CHANGE NAME OF SAVING FILE
#              # NB: CHECK X-AXIS VALUES --> MakeTimeSpan & plt.xlabel
#
#              stdev = (np.std(self.CalData[:,nChnl]))
##             str = "Standard deviation = %6.3f dB re [V]" % stdev
#              str = "$\sigma$ = %6.4f %s" % (stdev, self.Unit[nChnl]);
#              plt.title(self.basename + ', \n' + str)
#              str = "%s -> [%s] " %(self.Name[nChnl], self.Unit[nChnl]);
#              plt.ylabel(str) ; plt.autoscale(tight = True)
#              plt.xlabel('time ->[s]') #, \n NB: data = cutted pulses!
#              plt.grid()
#              plt.grid('on', 'minor')
#          plt.subplot(nPlot, 1, nSubplot-1)
##          plt.xlabel('time ->[s]'); # for common x-axis
#          plt.tight_layout()
#          plt.savefig(self.basename+'_data_1sec.png') #_data_1wav_zoom.png
#          plt.show()
#          plt.clf()
          
          
              ##################
              #### 0.1 sec ####
              ##################
              plt.plot(self.time[:self.Fs/10], self.CalData[:self.Fs/10,nChnl]) # <-- 1 window plotten
              plt.xlabel('#samples, \n NB: 1 sample = 5.12e-6 sec') # <<- x-axis for plotting 1 window
              # NB: CHANGE NAME OF SAVING FILE
              # NB: CHECK X-AXIS VALUES --> MakeTimeSpan & plt.xlabel

              stdev = (np.std(self.CalData[:,nChnl]))
#             str = "Standard deviation = %6.3f dB re [V]" % stdev
              str = "$\sigma$ = %6.4f %s" % (stdev, self.Unit[nChnl]);
              plt.title(self.basename + ', \n' + str)
              str = "%s -> [%s] " %(self.Name[nChnl], self.Unit[nChnl]);
              plt.ylabel(str) ; plt.autoscale(tight = True)
              plt.xlabel('time ->[s]')
              #plt.xlabel('time ->[s], nb: cutted!') #, \n NB: data = cutted pulses!
              plt.grid()
              plt.grid('on', 'minor')
          plt.subplot(nPlot, 1, nSubplot-1)
          #plt.xlabel('time ->[s]'); # for common x-axis
          #plt.tight_layout()
          plt.savefig(self.basename+'_data_0.1sec.png') #_data_1wav_zoom.png
          plt.show()

      else:
          nChnl = 0
          plt.plot(self.time, self.CalData[:, rChnl])
          stdev = (np.std(self.CalData[:, nChnl]))
#         str = "Standard deviation = %6.3f dB re [V]" % stdev
          str = "$\sigma$ = %6.4f %s" % (stdev, self.Unit[nChnl]);
          plt.title(str)
          str = "%s -> [%s] " %(self.Name[nChnl], self.Unit[nChnl]);
          plt.ylabel(str) ; plt.autoscale(tight = True)
          plt.grid()
          plt.grid('on', 'minor')
          plt.savefig(self.basename+'_data.png')
          plt.xlabel('time ->[s]');
          #plt.show()

     # plt.show()
  #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    def MakeSig(self, F0, nChnl):
        self.MakeTimeSpan()
        self.F0           = F0
        self.Data[:,nChnl]    = np.sin(2*pi*self.F0*self.time)

    def MakeNoise(self, std, nChnl):
        self.Data[:, nChnl]  = std*np.random.randn(self.length)

    def Spectrum(self, vecChnl):
        self.Calibrate()

        if self.PSDOn == 1: self.NFFT = self.Fs
        self.PSD        = np.zeros((self.NFFT/2,self.Nchnl), dtype=float)
        self.MakeWindow()
        self.MakeFreqSpan()
        self.length     = np.size(self.CalData,0)
        self.Segment    = int(self.length/self.NFFT)
        wnd             = self.window
        self.norm       = norm(wnd)**2
        double2single   = 2.0
        for nChnl in vecChnl:
            for span in range(0, self.Segment):
                Bg                 = span*self.NFFT
                end                = Bg+self.NFFT
                yw                 = wnd*self.CalData[Bg:end, np.int(nChnl)]
                a                  = fft(yw, self.NFFT)
                ac                 = np.conj(a)
                pxx                = np.abs(a*ac)
                self.PSD[:, nChnl] +=  double2single*pxx[0:self.NFFT/2]
#            self.Segment        = int(span/self.NFFT)+1
            self.PSD[:, nChnl]  /= (float(self.Segment)*self.NFFT*self.norm)
        return self.PSD
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    def TransferFunction(self, *args):
        self.Calibrate()
        nargin          = len(args)
        self.refr       = 0
        self.resp       = 1
        if nargin >= 1:
            self.refr = args[0]
        if nargin >= 2:
            self.resp = args[1]
        if self.PSDOn == 1: self.NFFT = self.Fs
        self.MakeWindow()
        self.MakeFreqSpan()

        self.length     = np.size(self.CalData, 0)
        self.Segment    = round(float(self.length)/float(self.NFFT))
        wnd             = self.window
        self.norm       = norm(wnd)**2
        Resp            = self.CalData[:, self.resp]
        Ref             = self.CalData[:, self.refr]
        nsamps          = self.NFFT
        start_index     = 0
        stop_index      = self.NFFT+start_index
        num_periods     = int(math.floor((len(Ref)-start_index)/self.NFFT))
        self.Sxx        = np.zeros((nsamps), dtype=float)
        self.Syy        = np.zeros((nsamps), dtype=float)
        self.Sxy        = np.zeros((nsamps), dtype=complex)
        for i in range(num_periods):
           respCalData  = wnd*Resp[i*self.NFFT+start_index:i*self.NFFT+stop_index]
           refCalData   = wnd*Ref[i*self.NFFT+start_index:i*self.NFFT+stop_index]
           self.Sx      = (2./self.NFFT)*np.fft.fft(refCalData)
           self.Sy      = (2./self.NFFT)*np.fft.fft(respCalData)
           self.Syy    += np.absolute(self.Sy)**2
           self.Sxx    += np.absolute(self.Sx)**2
           self.Sxy    += self.Sy*np.conjugate(self.Sx)

        WinNorm     = num_periods*self.norm    #Normalizing scale factor for window type
    #average the CalData
        self.Sxx   /= WinNorm
        self.Syy   /= WinNorm
        self.Sxy   /= WinNorm
        self.Txy    = self.Sxy/(self.Sxx+1e-15)
        self.Cxy    = np.absolute((self.Sxy)**2/(self.Sxx*self.Syy))
        self.Cxy[0] = 0 # DC -component has no useful info
        return self.Txy, self.Cxy
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    def PlotBode(self):

        Lrng     = int(len(self.Txy)/2)
        rng      = np.arange(Lrng, dtype=int)
        freqs    = np.fft.fftfreq(self.Txy.shape[0], d=(1./self.Fs))
        LowLim   = (float(self.Fs)/float(self.NFFT))
        UpperLim = max(freqs)
        plt.subplot(311)
        plt.semilogx(freqs[rng],20.*np.log10(np.absolute(self.Txy[rng])))
        plt.minorticks_on()
        plt.grid('on', which='both', axis='x')
        plt.grid('on', which='major', axis='y')
        plt.xlim(LowLim, UpperLim)
        plt.ylabel('Txy -> dB re []')

        plt.subplot(312)
        plt.semilogx(freqs[rng], 180.*((np.angle(self.Txy[rng])))/np.pi)
        plt.grid('on', which='both', axis='x')
        plt.grid('on', which='major', axis='y')
        plt.xlim(LowLim, UpperLim)
        plt.ylabel('Phase ->  [grd]')

        plt.subplot(313)
        plt.semilogx(freqs[rng], self.Cxy[rng])
        plt.grid('on',which='both', axis='x')
        plt.grid('on',which='major', axis='y')
        plt.xlim(LowLim, UpperLim)
        plt.ylabel('Coherence ->  []')
        plt.xlabel('Frequency -> [Hz]')
        plt.show()
 #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    def plotSpectrum(self, rChnl):
      nPlot    = len(rChnl)
      nSubplot = 1
      dF       = float(self.Fs)/float(self.NFFT)
      idxB     = int(self.FcumBegin/dF)
      idxE     = int(self.FcumEnd/dF)
      if self.OnePlot == 0:
          for nChnl in rChnl:
              plt.subplot(nPlot, 1, nSubplot)
              nSubplot += 1
              cummulative = self.PSD[idxB:idxE,nChnl].cumsum();
              std = 10*np.log10(cummulative[-1]); #According to Parseval
              plt.plot(self.freq, 10.*np.log10(abs(self.PSD[:, nChnl]))) # PDS in dB re ADC^2
              #plt.plot(self.freq, np.sqrt(self.PSD[:, nChnl])) #PSD in ADC^2
              # PSD -> abs -> ASD
              if self.Cummulative == 1:
                  plt.plot(self.freq[idxB:idxE],10.*np.log10(abs(cummulative)))
                  plt.plot(self.freq[idxE:idxB-1:-1], 10.*np.log10(abs(self.PSD[idxE:idxB-1:-1, nChnl].cumsum())))
              str = "$\sigma$ = %6.4f dB re %s" % (std,self.Unit[nChnl]);
              plt.title(self.basename + ', \n' + str)
              str = "%s -> [dB re %s] " %(self.Name[nChnl], self.Unit[nChnl]);
              if self.PSDOn == 1:str = "%s -> [dB re %s/Hz]" % (self.Name[nChnl], self.Unit[nChnl]);
              plt.ylabel(str) ; plt.autoscale(tight = True)
              plt.grid()
              plt.grid('on', 'minor')
          plt.subplot(nPlot,1, nSubplot-1)
          #plt.tight_layout()
          plt.xlabel('Frequency ->[Hz]\n ');
          plt.savefig(self.basename+'_spectrum.png')
          plt.show()
      else:
          plt.plot(self.freq, np.sqrt(abs(self.PSD[:, rChnl])))
          nChnl  = rChnl[0]
          str = "%s -> [%s] " %(self.Name[nChnl], self.Unit[nChnl]);
          if self.PSDOn == 1:str = "%s ->  [%s/$\sqrt{Hz}$]" % (self.Name[nChnl], self.Unit[nChnl]);
          plt.ylabel(str) ; plt.autoscale(tight = True)
#          plt.legend(self.Name[0])
#          plt.legend(self.Name[7])
          plt.grid()
          plt.grid('on','minor')
          plt.xlabel('Frequency ->[Hz]');
          #plt.show()

 #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    def plotSpectrumLinY(self, rChnl):
      nPlot    = len(rChnl)
      nSubplot = 1
      dF       = float(self.Fs)/float(self.NFFT)
      idxB     = int(self.FcumBegin/dF)
      idxE     = int(self.FcumEnd/dF)
      if self.OnePlot == 0:
          for nChnl in (rChnl):
              plt.subplot(nPlot, 1, nSubplot)
              nSubplot += 1
              cummulative = self.PSD[idxB:idxE, nChnl].cumsum();
              std = np.sqrt(cummulative[-1]); #According to Parseval
              plt.plot(self.freq, np.sqrt(abs(self.PSD[:, nChnl])))
              if self.Cummulative == 1:
                  plt.plot(self.freq[idxB:idxE], np.sqrt(abs(cummulative)))
                  plt.plot(self.freq[idxE:idxB-1:-1], np.sqrt(abs(self.PSD[idxE:idxB-1:-1, nChnl].cumsum())))
              str = "$\sigma$ = %6.4f %s" % (std, self.Unit[nChnl]);
              plt.title(str)
#        plt.title(r'$\alpha > \beta$')

              str = "%s -> [%s] " %(self.Name[nChnl], self.Unit[nChnl]);
              if self.PSDOn == 1:str = "%s ->  [%s/$\sqrt{Hz}$]" % (self.Name[nChnl], self.Unit[nChnl]);
              plt.ylabel(str) ; plt.autoscale(tight = True)
              plt.grid()
              plt.grid('on','minor')
              plt.show()

          plt.subplot(nPlot,1, nSubplot-1)
          plt.xlabel('Frequency ->[Hz]');
          plt.show()
      else:
          plt.plot(self.freq,np.sqrt(abs(self.PSD[:, rChnl])))
          nChnl  = rChnl[0]
          str = "%s -> [%s] " %(self.Name[nChnl], self.Unit[nChnl]);
          if self.PSDOn == 1:str = "%s ->  [%s/$\sqrt{Hz}$]" % (self.Name[nChnl], self.Unit[nChnl]);
          plt.ylabel(str) ; plt.autoscale(tight = True)
#          plt.legend(self.Name[0])
#          plt.legend(self.Name[7])
          plt.grid()
          plt.grid('on', 'minor')
          plt.xlabel('Frequency ->[Hz]');
          plt.show()

 #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    def plotSpectrumSingleLinY(self, nChnl):
      nSubplot = 1
      dF       = float(self.Fs)/float(self.NFFT)
      idxB     = int(self.FcumBegin/dF)
      idxE     = int(self.FcumEnd/dF)
      if self.OnePlot == 0:

          nSubplot += 1
          cummulative = self.PSD[idxB:idxE, nChnl].cumsum();
          std = np.sqrt(cummulative[-1]); #According to Parseval
          plt.plot(self.freq, np.sqrt(abs(self.PSD[:, nChnl])))
          if self.Cummulative == 1:
              plt.plot(self.freq[idxB:idxE], np.sqrt(abs(cummulative)))
              plt.plot(self.freq[idxE:idxB-1:-1], np.sqrt(abs(self.PSD[idxE:idxB-1:-1, nChnl].cumsum())))
          str = "$\sigma$ = %6.4f %s" % (std,self.Unit[nChnl]);
          plt.title(str)
#        plt.title(r'$\alpha > \beta$')

          str = "%s -> [%s] " %(self.Name[nChnl], self.Unit[nChnl]);
          if self.PSDOn == 1:str = "%s ->  [%s/$\sqrt{Hz}$]" % (self.Name[nChnl], self.Unit[nChnl]);
          plt.ylabel(str) ; plt.autoscale(tight = True)
          plt.grid()
          plt.grid('on','minor')
          plt.show()

#          plt.subplot(nPlot,1,nSubplot-1)
          plt.xlabel('Frequency ->[Hz]');
          plt.show()
      else:
          plt.plot(self.freq, np.sqrt(abs(self.PSD[:, nChnl])))
#          nChnl  = nChnl[0]
          str = "%s -> [%s] " %(self.Name[nChnl],self.Unit[nChnl]);
          if self.PSDOn == 1:str = "%s ->  [%s/$\sqrt{Hz}$]" % (self.Name[nChnl], self.Unit[nChnl]);
          plt.ylabel(str) ; plt.autoscale(tight = True)
#          plt.legend(self.Name[0])
#          plt.legend(self.Name[7])
          plt.grid()
          plt.grid('on', 'minor')
          plt.xlabel('Frequency ->[Hz]'); 
          plt.show()

 #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    def plotSpectrumLogLog(self, rChnl):
      nPlot    = len(rChnl)
      nSubplot = 1
      dF       = float(self.Fs)/float(self.NFFT)
      idxB     = int(self.FcumBegin/dF)
      idxE     = int(self.FcumEnd/dF)
      if self.OnePlot == 0:
          for nChnl in rChnl:
              plt.subplot(nPlot,1, nSubplot)
              nSubplot += 1
              cummulative = self.PSD[idxB:idxE, nChnl].cumsum();
              std = np.sqrt(cummulative[-1]); #According to Parseval
              plt.loglog(self.freq, np.sqrt(abs(self.PSD[:, nChnl])))
              if self.Cummulative == 1:
                  plt.loglog(self.freq[idxB:idxE], np.sqrt(abs(cummulative)))
                  plt.loglog(self.freq[idxE:idxB-1:-1], np.sqrt(abs(self.PSD[idxE:idxB-1:-1, nChnl].cumsum())))
              str = "$\sigma$ = %6.4f %s" % (std, self.Unit[nChnl]);
              plt.title(str)
#        plt.title(r'$\alpha > \beta$')

              str = "%s -> [%s] " %(self.Name[nChnl], self.Unit[nChnl]);
              if self.PSDOn == 1:str = "%s ->  [%s/$\sqrt{Hz}$]" % (self.Name[nChnl], self.Unit[nChnl]);
              plt.ylabel(str) ; plt.autoscale(tight = True)
              plt.grid()
              plt.grid('on','minor')
              plt.show()

          plt.subplot(nPlot, 1, nSubplot-1)
          plt.xlabel('Frequency ->[Hz]');
          plt.show()
      else:
          plt.loglog(self.freq, np.sqrt(abs(self.PSD[:, rChnl])))
          nChnl  = rChnl[0]
          str = "%s -> [%s] " %(self.Name[nChnl], self.Unit[nChnl]);
          if self.PSDOn == 1:str = "%s ->  [%s/$\sqrt{Hz}$]" % (self.Name[nChnl], self.Unit[nChnl]);
          plt.ylabel(str) ; plt.autoscale(tight = True)
#          plt.legend(self.Name[0])
#          plt.legend(self.Name[7])
          plt.grid()
          plt.grid('on', 'minor')
          plt.xlabel('Frequency ->[Hz]');
          plt.show()

 #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

    def plotSpectrumLog(self, rChnl):
      nPlot    = len(rChnl)
      nSubplot = 1
      dF       = float(self.Fs)/float(self.NFFT)
      idxB     = int(self.FcumBegin/dF)
      idxE     = int(self.FcumEnd/dF)

      for nChnl in rChnl:
        plt.subplot(nPlot,1, nSubplot)
        nSubplot += 1

        cummulative = self.PSD[idxB:idxE, nChnl].cumsum();
        std = 10*np.log10(cummulative[-1]); #According to Parseval
        plt.semilogx(self.freq, 10.*np.log10(abs(self.PSD[:,nChnl])))
        if self.Cummulative == 1:
            plt.semilogx(self.freq[idxB:idxE], 10.*np.log10(abs(cummulative)))
            plt.semilogx(self.freq[idxE:idxB-1:-1], 10.*np.log10(abs(self.PSD[idxE:idxB-1:-1, nChnl].cumsum())))
        str = "$\sigma$ = %6.4f dB re %s" % (std,self.Unit[nChnl]);
        plt.title(str)
        str = "%s -> [dB re %s] " %(self.Name[nChnl],self.Unit[nChnl]);
        if self.PSDOn == 1:str = "%s -> [dB re %s/$\sqrt{Hz}$]" % (self.Name[nChnl], self.Unit[nChnl]);
        plt.ylabel(str) ; plt.autoscale(tight = True)
        plt.grid()
        plt.grid('on','minor')
        plt.show()
      plt.subplot(nPlot, 1, nSubplot-1)
      plt.xlabel('Frequency ->[Hz]');
      plt.show()
  #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    def SetCumSumRange(self, Fbegin, Fend):
       if Fbegin > Fend:
           tmp    = Fbegin
           Fbegin = Fend
           Fend   = tmp
   
       if Fend > (self.NFFT/2):Fend = self.NFFT/2
       if Fbegin < 0:Fbegin = 0
       self.FcumBegin = Fbegin
       self.FcumEnd   = Fend
   #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

    def MakeWindow(self):
        WnTp   = self.WindowType.lower()
        if WnTp == 'rect':
            self.window = np.ones(self.NFFT)
        elif WnTp == 'hanning':
            self.window   = np.hanning(self.NFFT)
        elif WnTp == 'blackmanharris':
            self.window   = self.SpecialBlackmanHarris7T()
        else:
            self.window   = np.hanning(self.NFFT)

        return self.window

    def SetWindowType(self, type):
        self.WindowType = type
        return self.WindowType

    def SpecialBlackmanHarris7T(self):
    # coeff. 4-term window 5-term window 6-term window 7-term window
    # A0 3.635819267707608e-001 3.232153788877343e-001 2.935578950102797e-001 2.712203605850388e-001
    # A1 4.891774371450171e-001 4.714921439576260e-001 4.519357723474506e-001 4.334446123274422e-001
    # A2 1.365995139786921e-001 1.755341299601972e-001 2.014164714263962e-001 2.180041228929303e-001
    # A3 1.064112210553003e-002 2.849699010614994e-002 4.792610922105837e-002 6.578534329560609e-002
    # A4 1.261357088292677e-003 5.026196426859393e-003 1.076186730534183e-002
    # A5 1.375555679558877e-004 7.700127105808265e-004
    # A6 1.368088305992921e-005

        A0 =2.712203605850388e-001
        A1 =4.334446123274422e-001
        A2 =2.180041228929303e-001
        A3 =6.578534329560609e-002
        A4 =1.076186730534183e-002
        A5 =7.700127105808265e-004
        A6 =1.368088305992921e-005
        N  = self.NFFT
        n  = 2.*np.linspace(0, N-1, N)/N
        self.window = A0-A1*cos(pi*n)+A2*cos(2*pi*n)-A3*cos(3*pi*n)+A4*cos(4*pi*n)-A5*cos(5*pi*n)+A6*cos(6*pi*n)
        return self.window

    def MakeFreqSpan(self):
        self.freq     = np.linspace(0, self.NFFT/2-1, self.NFFT/2)*float(self.Fs)/float(self.NFFT)
        return self.freq

    def MakeTimeSpan(self):
        self.time     = np.linspace(0, self.length-1, self.length)/self.Fs
        return self.time

    def Decimate(self, Factor):
#        FNyq        = 0.5*self.Fs/Factor       #New Nyquist Freq
#        F0          = 0.95*FNyq               #New Bandpass Freq
        tmp         = signal.decimate(self.Data, Factor, n=None, ftype='fir', axis=0)
        self.Data   = tmp
        self.length = len(tmp)
        self.Fs    /= Factor

    def ShiftLeftSample(self, rChnl, No):

        for nChnl in rChnl:
            self.Data[:, nChnl] = np.roll(self.Data[:, nChnl], No)
