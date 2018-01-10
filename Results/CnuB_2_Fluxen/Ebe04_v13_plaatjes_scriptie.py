from math import *
import matplotlib.pyplot as plt
import numpy as np
import sys
import scipy.integrate as integrate
import ConfigParser
#from Statistics_config import *

# ============================================================================
# MAIN ARTICLE
# ============================================================================
# Ebe04: 
#   Eberle 2004
#   Relic Neutrino Absorption Spectroscopy
#   arxiv:hep-ph/0401203
# ============================================================================
# OTHER ARTICLES referred to in comments:
# ============================================================================
# Aar15:
#   Aartsen, 2015
#   Atmospheric and Astrophysical Neutrinos above 1 TeV Interacting in IceCube
#   arXiv:1410.1749v2
# Ahl14:
#   Ahlers, 2014
#   Pinpointing extragalactic neutrino sources in light of recent IceCube 
#       observations
#   arXiv: 1406.2160v1
# Bau17:
#   Baumann, 2017
#   Syllabus Cosmology course UvA 2017
#   http://www.damtp.cam.ac.uk/user/db275/Cosmology/Lectures.pdf
# Fan16:
#   Fang, 2016
#   High-energy neutrinos from sources in clusters of galaxies
#   arXiv:1607.00380
# Fod03:
#   Fodor, 2003
#   Bounds on the cosmogenic neutrino flux
#   arXiv:hep-ph/0309171v1
# Oli06:
#   D'Olivo, 2006
#   UHE neutrino damping in a thermal gas of relic neutrinos
#   arXiv:astro-ph/0507333v2
# Lun13:
#   Lunardini 2013
#   Ultra High Energy Neutrinos: Absorption, Thermal Effects and Signatures
#   arXiv:1306.1808v2
# Wic04:
#   Wick, 2004
#   High-energy cosmic rays from gamma-ray bursts
#   arXiv:astro-ph/0310667v2
# ============================================================================


# =================
# Parameter configuration
# =================

def Config(var):
    config = ConfigParser.RawConfigParser()
    config.read('Ebe04_Flux_setup.cfg')
    
    if var == 'm_n': return config.getfloat('HE_Neutrinos', 'm_n')
    elif var == 'p_n': return config.getfloat('HE_Neutrinos', 'p_n')
    elif var == 'num_dens': return config.getfloat('Relic_Neutrinos', 'num_dens')
    elif var == 'm_Z': return config.getfloat('Z_boson', 'm_Z')
    elif var == 'gamma_Z': return config.getfloat('Z_boson', 'gamma_Z')
    elif var == 'h': return config.getfloat('Universe', 'h')
    elif var == 'Om': return config.getfloat('Universe', 'Om')
    elif var == 'Ol': return config.getfloat('Universe', 'Ol')
    elif var == 'Ok': return config.getfloat('Universe', 'Ok')
    elif var == 'z_min': return config.get('Universe', 'z_min')
    elif var == 'z_max': return config.get('Universe', 'z_max')
    elif var == 'n': return config.get('Universe', 'n')
    elif var == 'alpha': return config.get('Universe', 'alpha')
    elif var == 'eta0': return config.getfloat('Source_Emissivity', 'eta0')
    elif var == 'j': return config.getfloat('Source_Emissivity', 'j')
    else: print(var, " is an unknown parameter"), sys.exit()


# =================
# "Basic formula's"
#  - Resonance energy
#  - Hubble parameter
# =================

def E_resonance(m_n):
    ''' Function to determine the resonance energy [eV] 
    for given neutrino mass [eV]'''
    m_z = 91.2e9    # [eV]
    return m_z**2 / (2. * m_n)                                                  # (eq 1)

def Hubble(z, h=0.678, Om=0.308, Ol=0.692, Ok=0.):
    '''Function to determine the hubble parameter at given z. The curvature 
    parameters are given a default values resulting in a flat universe'''       # default values from Ebe04
    H0 = 100 * h    # 100 h [km/s/Mpc]                                          <== NOG KEER EEN CONSTANTE VOOR DIMENSIONS, however, J and eta0 in source emisivity also unknown --> doesn't matter
    return sqrt(H0**2 * ( Om * (1+z)**3 + Ok * (1+z)**2 + Ol) )                 # (eq 8)


# =================
# Crossection & Annihilation Probabiity
# =================

def Oli05_crossection_Sigma(e):                                                 # NB: ultra relativistic approximation 
    '''Funtion that returns the crossection in [nb], given the UHE neutrino 
    energy 'e' in [eV]'''
    # D'Olivo et al 2005, equation 22
    m       = 0.1               # neutrino mass                 in [eV]
    P       = 6.08e-4           # neutrino momenta              in [eV]
    GammaZ  = 2.4952*1e9        # width for Z-decay to fermions in [eV]
    MassZ   = 91.1876*1e9       # Z-mass                        in [eV]
    GF      = 1.166364e-23      # Fermi coupling constant       in [eV^-2]

    Ep      = sqrt(P**2 + m**2)
    xi = GammaZ**2/MassZ**2
    teller = (1+xi)*4*e**2*(Ep+P)**2 - 4*MassZ**2*e*(Ep+P) + MassZ**4
    noemer = (1+xi)*4*e**2*(Ep-P)**2 - 4*MassZ**2*e*(Ep-P) + MassZ**4

    aid1 = (2*e*(1+xi)*(Ep+P)-MassZ**2)/(GammaZ*MassZ)
    aid2 = (2*e*(1+xi)*(Ep-P)-MassZ**2)/(GammaZ*MassZ)

    # convert cross section units:
    C = (1/(2.56819e-24))       # 1 [nb] = 1e-37 [m^2] = 2.56819e-24 [eV^-2]

    return ( C * (2*sqrt(2.)*GF * GammaZ * MassZ) / (e * Ep) * \
             (
                 (1/(1 + xi))
                 + MassZ**2/(4*e*P*(1+xi)**2) * log(teller/noemer) \
                 + ((1-xi)*MassZ**3 / ((1+xi)**2 * 4*e*P*GammaZ)) * \
                 ( atan(aid1) - atan(aid2))
              )
            )

def Ebe04_Oli05_annihilation_probability(e, h):
    num_dens = 56 # in [cm^-3]
    sigma = Oli05_crossection_Sigma(e)
    H0 = Hubble(0)
    
    # [cm^-3 * nb * km^-1 * s * Mpc] to seconds and meters
    C1 = 1e-37 * 3.085e22 / (1e-6 * 1e3)
    # [s]/[m] to natual units
    C2 = 0.197e-15/6.58e-25

    fraction = num_dens * sigma / H0 * C1 * C2

    return 0.71/h * fraction

# =================
# Survival Probability
# =================

# -- as a function of E/Eres and z
def Ebe04_survival_probability_P(e, Eres, z, 
                                 h=0.678, Om=0.308, Ol=0.692, Ok=0.):
    '''Funtion that returns survival probability for the e values 
    within range: Eres/(1+z) < e < Eres ''' 

    ann_prob = 0.71/h * 0.03        # anihilation probability                   # (Ebe04: eq16)
#    ann_prob = Ebe04_Oli05_annihilation_probability(e, h)                      # Oli06: annihilation probability (ann_prob) better if crossections as a function of energy
    assert(float(Om + Ol + Ok) == 1.), "No Flat Universe!"

    if e < Eres/(1+z) or e > Eres:
        return 1.                   # No absorption
    else:
        return exp(-ann_prob * (((Eres/e)**3) / sqrt(Om*((Eres/e)**3) 
        + Ok*((Eres/e)**2) + Ol)))                                              # (eq15)


# =================
# source emissivity distribution
#   NB: the function has been split up because there is a bug when integrating 
#   over the non-z-dependent factors.
# =================

def Ebe04_source_emissivity_L_z(z, z_max, n_min_alpha):                         # Wic04: explaines powerlaw ansatz. NB: approximation is appropriate upto z ~ 2, but for E_{CR} < 10^18 eV -->  predicting too high flux
    '''Function that returns all parts of the source emissivity function that
    DO depend on z'''
    z_min = 0
    if z > z_min:
        if z_max > z:
            return (1 + z) ** n_min_alpha                                       # (eq 24 --> eq 27&28, but only z-dependent part (for integration))
        else: return 0.
    return 0.

def Ebe04_source_emissivity_L_no_z(e, alpha):                                   # VERY UNCERTAIN!
    '''Function that returns all parts of the source emissivity function that
    DO NOT depend on z'''
    eta0 = 1    # 1e-5 [Mpc^-3]                                                 # Ahl14: local source density in [Mpc^-3] different units...
    j = 1                                                                       # normalisation factor ?
    return eta0 * j * e**(-alpha)                                               # (eq 24 --> eq 27&28, but only z-independent part)


# =================
# neutrino flux
# - primary & secondary
# =================

def Ebe04_neutrino_flux_earth_F(e, Eres, z_max, n, alpha, absorption=True):
    '''Function to determine the Flux per neutrino flavor by integrating 
        the survival probability * source emissivity / Hubble parameter
        over all redshifts z and multiply it by some constants.
        When no absorption is included the survival probabiliy is set to 1'''

    if absorption == True:
        # - function to integrate for Flux
        f1_to_integrate = lambda z: 1/Hubble(z) * \
                            Ebe04_survival_probability_P(e, Eres, z) * \
                            Ebe04_source_emissivity_L_z(z, z_max, n-alpha)
    elif absorption == False:
        # - function to integrate for Flux without absorption
        f1_to_integrate = lambda z: 1/Hubble(z) * \
                            Ebe04_source_emissivity_L_z(z, z_max, n-alpha)
                                                                                # INTEGRATION BOUNDARIES:
    integrant1, err1 = integrate.quad(f1_to_integrate, 0, z_max)                # if z > z_max; L (Ebe04_source_emissivity_L_z) returns 0, therfore this doesn't add anything to the integral
                                                                                # Lun13: for E > 10^11 GeV, neutrinio horizon of z ~ 140, beyond which universe is opaque to neutrinos <= taking this as upperlimit gives the same
                                                                                # Bau17: neutrino decoupling at redshift z = 6e9 <= taking this as upperlimit goes wrong... ?
                                                                                # Ebe04: in eq 23 an integral to np.inf is shown. But this gives the same result as z=z_max, and z=140 since z>z_max doesn't add anything
    return 1/(4 * pi) * integrant1 * 1/3. * Ebe04_source_emissivity_L_no_z(e, 
            alpha)


def Ebe05_neutrino_flux_earth_F(e, Eres, z_max, n, alpha, Z_decay=False, \
                                absorption=True):
    '''Function to determine the Flux per neutrino flavor with the option to
    include primary flux only (Z_decay = False),
    or include primary flux as well as secondary flux (Z_decay = True)'''

    primary_flux = Ebe04_neutrino_flux_earth_F(e, Eres, z_max, n, alpha, 
                                               absorption)

    if absorption == False:
        # if no absorption; no UHE Z-boson created, so no secondary neutrinos
        return primary_flux

    if Z_decay == False:
        return primary_flux

    elif Z_decay == True:
        decay_frac_to_nu = 0.4                                                  # Z -> nu nu in 20% of the decaying processes. here times 2 since doubling of neutrino detection possibility
        secondary_flux = Ebe04_neutrino_flux_earth_F(2*e, Eres, z_max, n, 
                                                  alpha, absorption=False) - \
                         Ebe04_neutrino_flux_earth_F(2*e, Eres, z_max, n, 
                                                  alpha, absorption=True)
        return primary_flux + secondary_flux * decay_frac_to_nu


# =================
# Plot functions
# =================

def plot_sigma():
    Eres = E_resonance(0.1)
    E = np.logspace(-2, 3, 10000)*Eres
    
    plt.plot(E, np.array([Oli05_crossection_Sigma(x) for x in E]), '-')
    plt.yscale('log')
    plt.xscale('log')
    plt.xlabel("E$\\nu$ [eV]", size =14, position=(1, 0), ha='right', fontsize=14)
    plt.ylabel("$\sigma$ [nb]", size =14, position=(0,1), ha='right', fontsize=14)
    plt.show()

def plot_P():
    redshift = [2., 5., 20.]
    Eres = E_resonance(0.1)
    for z in redshift:
        if z == 2: col = 'r'
        elif z == 5: col = 'b'
        elif z == 20: col = 'g'
        E = np.linspace(Eres/(1-z), Eres*1.5, 1000)
        P = []
        for e in E:
            P.append(Ebe04_survival_probability_P(e, Eres, z))
        labeltxt = "$\mathrm{z} = %i$" %z
#        plt.plot(E, P, "-", ms = 5, color = col, label = labeltxt)
        plt.plot(E/Eres, P, "-", ms = 5, color = col, label = labeltxt)
    plt.legend(loc = 'lower left')
    plt.xscale('log')
#    plt.xlabel("E (eV)", fontsize=14)
    plt.xlabel("E/E$^{\mathrm{res}}$", fontsize=14)
    plt.ylabel("P (E(1+z), z)", fontsize=14)
#    plt.xlim(0.01*Eres, 1.5*Eres)
    plt.xlim(0.01, 1.5)
    plt.ylim(0, 1.1)
    plt.title('Survival Probability',fontsize=18) #, Ebe04
    plt.text(2e22, 0.25, r'$ m_{\nu} = 0.1 eV $', size = 18)
    plt.show()


def plot_F(frac = False, Z_decay = False):
    '''Function to plot the flux.
    If frac = True, the flux will be plotted divided by the flux when there
    is no absorption.'''
                                                                                # Fod03, fitting procedure gives the most probable values for Emax, alpha and n.
    redshifts_max = [2, 20]                                           # z_max
    alphas        = [2]                                       # spectral index -> (Ebe04: [1-2])
    ns            = [2, 4]                                              # powers of (1+z) in activity Ebe04 eq27

    for z_max in redshifts_max:
        if z_max == 2: col = 'r'
        elif z_max == 20: col = 'g'

        for n in ns:
            for alpha in alphas:
                if n - alpha == 2 and z_max == 2:
                    labeltxt = 'Bottom up'
                    linst = '-'
                elif n - alpha == 0 and z_max == 20:
                    labeltxt = 'Top Down'
                    linst = '-'
                else: continue # <= for now only analyse two scenario's

                Eres = E_resonance(0.1)
                E = np.logspace(-2, 3, 10000)*Eres                              # from 10^a to 10^b in c logarithmic steps
                F = []
                F1 = []

                if frac == True:

                    Fnoabs = []
                    F_frac_Fnoabs = []
                    F1noabs = []
                    F1_frac_Fnoabs = []

                    for e in E:
                        F.append(Ebe05_neutrino_flux_earth_F(e, Eres, z_max, n, 
                                      alpha, Z_decay, absorption = True))
                        Fnoabs.append(Ebe05_neutrino_flux_earth_F(e, Eres, 
                                      z_max, n, alpha, Z_decay, 
                                      absorption = False))
                        F_frac_Fnoabs.append(F[-1]/Fnoabs[-1])

#                        if Z_decay == True:
#                        # If flux + secondary flux is plotted: Z_decay = True,
#                        # add plot of primary flux only in black, to compare:
#                            F1.append(Ebe05_neutrino_flux_earth_F(e, Eres, 
#                                      z_max, n, alpha, Z_decay=False, 
#                                      absorption = True))
#                            F1noabs.append(Ebe05_neutrino_flux_earth_F(e, Eres, 
#                                      z_max, n, alpha, Z_decay=False,
#                                      absorption = False))
#                            F1_frac_Fnoabs.append(F1[-1]/F1noabs[-1])

                    labeltxt = r"$\mathrm{z_{max}} = %i,\ \mathrm{n} - \alpha = %i$" %(z_max, n-alpha)
                    print 'plot', labeltxt
                    plt.plot(E/Eres, F_frac_Fnoabs, linestyle = linst, 
                             color = col, label = labeltxt)
#                    if Z_decay == True: 
#                        plt.plot(E/Eres, F1_frac_Fnoabs, linestyle = linst, 
#                             color = 'black', label = labeltxt)
                    plt.ylabel("F/F$_{\mathrm{no\ abs}}$", fontsize=14)
                    plt.text(5e20, 0.65, r'$ m_{\nu} = 0.1 eV $', size = 18)
                    plt.ylim(0, 1.1)
                    plt.xlim(1e-2, 1.5)
                    plt.xlabel("E/E$^{\mathrm{res}}$", fontsize=14)


                else:

                    for e in E:
                        F.append(Ebe05_neutrino_flux_earth_F(e, Eres, z_max, n, 
                                     alpha, Z_decay, absorption = True) * e**2)       # <======================================  F * E^2
                    print 'plot', labeltxt
                    print "n = %i, alpha = %i, z_max = %i" %(n, alpha, z_max)
                    #labeltxt = "z_max = %i, n = %i, alpha = %i" %(z_max, n, alpha)
                    plt.plot(E, F, linestyle = linst, color = col, 
                             label = labeltxt)
 #                   plt.ylabel("F(E) [Arbitrary Units]", fontsize=14)                  # <====================================== LABEL
                    plt.ylabel(r"F(E) $\cdot$ E$^{2}$ [Arbitrary Units]", fontsize=14)  # <====================================== LABEL
                    plt.yscale('log')                                                   # <====================================== LOG
                    #plt.text(5e21, 1.7e-46, r'$ m_{\nu} = 0.1 eV $', size = 18)
                    plt.xlim(4e20, 4e23)
                    plt.xlabel("E [eV]", fontsize=14)

    plt.legend(loc = 'best')
    if Z_decay == True:
#        plt.title("Normalized Flux with Z-decay", fontsize=18)#, Ebe04+Ebe05")
        plt.title("Neutrino Flux at Earth", fontsize=18)#, Ebe04+Ebe05")
    elif Z_decay == False:
        plt.title("Normalized Flux", fontsize=18)#, Ebe04")
    plt.xscale('log')                                                                  # <======================================  LOG
    plt.show()


# =================
# write data to txt file
# =================

def flux_in_txt(m_n , z_max, n, alpha, Z_decay = False, cros = 'var'):

    Eres = E_resonance(m_n)
    #E = np.logspace(-3, 3, 10000)*Eres
    E = np.linspace(1e18, 5e22, 10000)

    H0savename = "FluxE2_H0data_m%s_zmax%i_n%i_alpha%i_Zdecay%s_cros%s.txt" %(m_n , z_max, n, alpha, Z_decay, cros)
    H1savename = "FluxE2_H1data_m%s_zmax%i_n%i_alpha%i_Zdecay%s_cros%s.txt" %(m_n , z_max, n, alpha, Z_decay, cros)

    ###########################################################################
    if H0savename in  os.listdir(os.curdir):
        print "Do you want to overwrite the files \n - %s \n - %s" %(H0savename, H1savename)
        var = raw_input("Please enter Y or N: ")
        while var != 'Y' and var != 'N':
            print "###\nUnknown input, please try again:"
            print "Do you want to overwrite the files \n - %s \n - %s" %(H0savename, H1savename)
            var = raw_input("Please enter Y or N: ")
        if var == 'N':
            sys.exit()
    ###########################################################################

    # H0 data: there is NO dip in the spectrum
    print "Writing data to txt file. One moment please..."
    with open(H0savename, 'w') as f:
        for e in E:
            FluxE2 = Ebe05_neutrino_flux_earth_F(e, Eres, z_max, n, alpha, Z_decay, absorption = False)*(e**2)
            line = "%s,%s\n" %(e, FluxE2)
            f.write(line)

    # H1 data: there is A dip in the spectrum
    with open(H1savename, 'w') as f:
        for e in E:
            FluxE2 = Ebe05_neutrino_flux_earth_F(e, Eres, z_max, n, alpha, Z_decay, absorption = True)*(e**2)
            line = "%s,%s\n" %(e, FluxE2)
            f.write(line)
# =================
# main
# =================

def main():

    # -- plot crossections
#    plot_sigma()

    # -- plot absorption probability as a function of Energy
#    plot_P()

    # -- plot Flux as a fraction of Flux in case of no absorbtion, 
    #    as a function of Energy/Eres
#    plot_F(frac = True, Z_decay = False)
#    plot_F(frac = True, Z_decay = True)
                                                                                # Aar15 referred to in Fan16:
    # -- plot Flux as a function of Energy                                      # The all-flavor diffuse neutrino flux is reported to be phi = 2.06e-18 (E / (1e5 GeV))**-2.46 GeV^-1 cm^-1 sr^-1 s^-1
#    plot_F(Z_decay = False) # F * E^2 is plotted here                          # for the energy range of 25TeV < E < 1.4 PeV (Aartsen et al. 2015)
    plot_F(Z_decay = True)  # F * E^2 is plotted here

#    flux_in_txt(m_n=0.1 , z_max=20, n=6, alpha=2, Z_decay = True, cros='Const')

if __name__ == "__main__":
    sys.exit(main())

