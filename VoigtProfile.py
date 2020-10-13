import numpy as np
from scipy.special import wofz
from scipy.constants import c, pi, m_e, elementary_charge
import astropy.constants as cst
import matplotlib.pyplot as plt

def VoigtProfile(x, sigma, gamma):
    """
    Function to return the value of a (normalized) Voigt profile centered at x=0
    and with (Gaussian) width sigma and Lorentz damping (=HWHM) gamma. 

    The Voigt profile is computed using the scipy.special.wofz, which returns
    the value of the Faddeeva function.

    WARNING
    scipy.special.wofz is not compaible with np.float128 type parameters.

    Args:
        x (float64): Scalar or array of x-values
        sigma (float64): Gaussian sigma component
        gamma (float64): Lorentzian gamma (=HWHM) component

    Returns:
        ndarray: Flux array for given input

    """
    
    z = (x + 1j * gamma)/sigma/np.sqrt(2)

    return (
        np.real(wofz(z)) / sigma / np.sqrt(2*np.pi)
    )



def VoigtProfile_nu(nu, nu0=0, b=0, gamma=0):
    """
    Function to return the value of a (normalized) Voigt profile, expressed per unit bandwidth in frequency  ($\phi(\nu)$) 
    centered at frequency nu0. b is the broadening parameter, and gamma the Lorentz damping (=HWHM).  

    This is a wrapper around VoigtProfile that produces a normalized Voigt profile at the origin. 
    This method just does all the conversions. 

    Args:
        nu (float64): Scalar or array of frequency values [Hz]
        nu0 (float64): rest frequency of transition [Hz]
        b (float64): Thermal & turbulence broadening parameter [km/s]
        gamma (float64): Lorentzian gamma (=HWHM) 

    Returns:
        ndarray: Flux array for given input

    """

    # 1e3 to convert b from km/s to m/s
    Delta_nud = nu0/c * (b * 1.e3)      
    sigma = Delta_nud / np.sqrt(2.)
    x = nu-nu0
    #print("b = ", b)
    #print("sigma = ", sigma)

    return (
        VoigtProfile(x,sigma,gamma/(4*pi))
    )


def VoigtProfile_lambda(wave, lambda0=0, b=0, gamma=0):
    """
    Function to return the value of a (normalized) Voigt profile, expressed per unit wavelength [in AA]  ($\phi(\lambda)$) 
    centered at wave0. b is the broadening parameter, and gamma the Lorentz damping (=HWHM).  

    This is a wrapper around VoigtProfile_nu; here we just convert from wavelength to frequency: 
    \phi(\lambda) = \phi(\nu)d\nu/\dlambda

    Args:
        wave (float64): Scalar or array of wavelength values [AA] 
        lambda0 (float64): rest wavelength of transition [AA]
        b (float64): Thermal & turbulence broadening parameter [km/s]
        gamma (float64): Lorentzian gamma (=HWHM) units??

    Returns:
        ndarray: Flux array for given input

    """
    
    # Since wave is in AA, need to be careful! 
    nu = 1e10 * c / wave
    nu0 = 1e10 * c / lambda0

    return (
        VoigtProfile_nu(nu,nu0=nu0,b=b,gamma=gamma) * c / (wave * 1e-10)**2
    )


def tau_lambda_Voigt(wave, lambda0=0, b=0, gamma=0, N=0, f=0):
    """
    Function to return the value of a Voigt-shaped optical depth profile, expressed per unit wavelength  ($\phi(\lambda)$) 
    centered at wave0. b is the broadening parameter, and gamma the Lorentz damping (=HWHM); N the column density and f
    the oscillator strength.   

    This is a wrapper around VoigtProfile_lambda -- we just have to multiply by the proper scale factor and take care
    of units. 

    Args:
        wave_lambda (float64): Scalar or array of wavelength values [AA]
        lambda0 (float64): rest wavelength of transition [AA]
        b (float64): Thermal & turbulence broadening parameter [km/s]
        gamma (float64): Lorentzian gamma (=HWHM) units??
        N: column density (cm^{-2})
        f: dimensionless

    Returns:
        ndarray: Flux array for given input

    """

    factor = (np.sqrt(pi) * elementary_charge**2)/(m_e * c * (b*1e3)) * f * (1e-10 * lambda0) * (N*1e4)
    print(pi,elementary_charge,m_e,c,f,lambda0,N)

    return (
        factor * VoigtProfile_lambda(wave,lambda0=lambda0,b=b,gamma=gamma)
    )

def integrate_under_profile(x,y,string=None):
    """
    This function compares three different ways to calculate the integral of the function contained in x and y
    and prints out the results. 
    """
    Area_trap = trapz(y,x)
    Area_simp = simps(y,x)
    Area_idl = idl_tabulate(x,y)
    print("Results for " + string + ": Area = ", Area_trap, " | " , Area_simp, " | ", Area_idl)


if __name__ == "__main__":

    import scipy
    from scipy.integrate import trapz
    from scipy.integrate import simps
    from math import *
    from idl_tabulate import idl_tabulate
    from random import uniform

    """
    This is a series of tests to see if the implementation of the various Voigt functions has been done correctly. 
    It includes tests for the basic Voigt function, as well as for the various forms of the normalized Voigt profiles. 
    """
    
    # First tests: verify that VoigtProfile is properly normalized, and functionally correct.  
    # We will use a random value for sigma, and sample 100 sigma away on each side of the peak to get 
    # a large enough integration area. 
    sigma = uniform(1,5)
    print("Sigma =", sigma)
    step = 0.01
    x=np.arange(-100.*sigma,+100.*sigma,step)
    
    # For gamma=0, we should get a normalized Gaussian back. 
    GaussProf = VoigtProfile(x,sigma,0.)
    # Peak for a Gaussian should be at 1 / sigma / sqrt(2\pi)
    GaussPeak = 1./sigma/np.sqrt(2.*np.pi)
    GaussMax = GaussProf.max()
    print("Test for gamma=0: Maximum should be at", GaussPeak, "and is at", GaussMax, "so delta =", GaussPeak-GaussMax)
    integrate_under_profile(x,GaussProf,"Gaussian")
    # Let's compare to Lorentz profile with same FWHM, and create Voigt profiles with inbetween values. 
    GaussFWHM = sigma * 2. * sqrt(2.*log(2.))
    LorentzProf = VoigtProfile(x,1.e-10,GaussFWHM/2)
    LorentzPeak = 2. / (pi * GaussFWHM)
    integrate_under_profile(x,LorentzProf,"Lorentz")
    # Note that Lorentzian has broad wings -- so we don't get 1 typically. 
    VoigtProf1 = VoigtProfile(x,sigma,GaussFWHM/10)
    integrate_under_profile(x,VoigtProf1,"Voigt Profile 1")
    VoigtProf2 = VoigtProfile(x,sigma/2,GaussFWHM/2)
    integrate_under_profile(x,VoigtProf2,"Voigt Profile 2")
    
    # Now let's plot things
    plt.xlim(-10*sigma,+10*sigma)
    plt.plot(x, GaussProf)
    plt.plot(x, LorentzProf, color='orange')
    plt.plot(x, VoigtProf1, color='green')
    plt.plot(x, VoigtProf2, color='red')
    plt.axvline(x=-GaussFWHM/2,linestyle='--')
    plt.axvline(x=+GaussFWHM/2,linestyle='--')
    plt.axhline(y=GaussPeak,linestyle='--')
    plt.axhline(y=GaussPeak/2,linestyle='--')
    plt.axhline(y=LorentzPeak,linestyle='--', color='orange')
    plt.axhline(y=LorentzPeak/2,linestyle='--', color='orange')
    plt.show()

    # Now on to VoigtProfile_nu and VoigtProfile_lambda
    # To properly sample things, create velocity array in steps of 0.05 km/s. 
    vel = np.arange(-10.,10.,.05)
    # Convert to wavelength now, and let's center at KI line -- 7698.974 AA.  
    wave0=7698.974
    f=3.393e-1
    N=44.3e10
    b=0.72
    gamma=3.8e7
    delta_wave = wave0 * (vel * 1e3) / c
    wave = wave0 + delta_wave
    # Let's not forget that wave is in AA!
    # Calculate frequency grid 
    nu0 = 1e10 * c / wave0
    nu = 1e10 * c / wave
    # Let's create a profile with a b of 0.72 km/s. 
    profile_nu = VoigtProfile_nu(nu,nu0=nu0,b=b,gamma=gamma)
    profile_lambda = VoigtProfile_lambda(wave,lambda0=wave0,b=b,gamma=gamma)
    # Frequency order is reversed, so integrate over -nu
    integrate_under_profile(-nu,profile_nu,"Voigt_nu")
    # Make sure to integrate in MKS units! Wavelength is in AA!! 
    integrate_under_profile(wave*1e-10,profile_lambda,"Voigt_lambda")
        

    print(profile_nu.min(), profile_nu.max())
    plt.xlim(-5,5)
    plt.xlabel('Velocity [km/s]')
    plt.ylabel("Phi(nu)")
    plt.title("Phi(nu) versus velocity")
    plt.plot(vel,profile_nu)
    plt.show()
    #plt.plot(nu,profile_nu)
    #plt.show()
    #plt.xlim(wave0-.1,wave0+.1)
    plt.xlabel("Wavelength (AA)")
    plt.ylabel("Phi(lambda)")
    plt.title("Phi(lambda) versus velocity")
    plt.plot(vel,profile_lambda)
    print(profile_lambda.min(), profile_lambda.max()*1e-9)
    plt.show()

    # Finally, Voigt optical depth profile. 
    # Let's use Dan Welty's values for the Na I D lines (mean of hyperfine components)
    tau = tau_lambda_Voigt(wave, lambda0=wave0, b=b, gamma=gamma, N=N, f=f)
    plt.xlim(wave0-1,wave0+1)
    plt.xlabel("Wavelength (AA)")
    plt.ylabel("Optical Depth")
    plt.plot(wave,tau)
    plt.show()
    integrate_under_profile(wave*1e-10,tau,"tau_lambda")
    #plt.plot(wave,np.exp(-tau))
    #plt.show()