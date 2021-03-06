import numpy as np
import matplotlib.pyplot as plt
from edibles import PYTHONDIR
from edibles.utils.edibles_oracle import EdiblesOracle
from edibles.utils.edibles_spectrum import EdiblesSpectrum
from edibles.utils.voigt_profile import *
from pathlib import Path
import astropy.constants as cst
from scipy.interpolate import interp1d
from scipy.stats import pearsonr


if __name__ == "__main__":
   #############################################################
   #
   # EXAMPLE for HD 183143
   #
   #############################################################

   # Get one of the spectra -- this is just the 3302 region. 
   pythia = EdiblesOracle()
   List = pythia.getFilteredObsList(object=["HD 183143"], MergedOnly=True, Wave=3302.0)
   test = List.values.tolist()
   filename = test[0]
   print(filename)
   wrange = [3301.5, 3304.0]
   sp = EdiblesSpectrum(filename)
   wave = sp.wave
   # What would be the velocity resolution? Assume oversampling rate of 2. 
   print("Velocity resolution:", 2. * (wave[1]-wave[0])/wave[0] * cst.c.to("km/s").value)
   flux = sp.flux
   idx = np.where((wave > wrange[0]) & (wave < wrange[1]))
   wave = wave[idx]
   flux = flux[idx]
   flux = flux / np.median(flux)
   
   # Now create a model for a single cloud, at v_rad = 0. Let's pick N and b such that they produce a decent absorption line. 
   lambda0 = [3302.369, 3302.978]
   f = [8.26e-03, 4.06e-03]
   gamma = [6.280e7, 6.280e7]
   b = [1.0]
   N = [1e13]
   v_rad = [0.0]
   v_resolution = 4.00

   AbsorptionLine = voigt_absorption_line(
    wave,
    lambda0=lambda0,
    b=b,
    N=N,
    f=f,
    gamma=gamma,
    v_rad=v_rad,
    v_resolution=v_resolution,
    )
   #plt.plot(wave, flux)
   #plt.gca().get_xaxis().get_major_formatter().set_useOffset(False)
   #plt.plot(wave, AbsorptionLine, color="orange", marker="*")
   #plt.show()
   v_rad_grid = np.arange(-10.,40.,.5) # in km / s
   all_corr = v_rad_grid * 0.
   #print(v_rad_grid)
   for loop in range(len(v_rad_grid)):
      v_rad = v_rad_grid[loop]
      Doppler_factor = 1. + v_rad / cst.c.to("km/s").value
      #print(Doppler_factor)
      new_wave = wave * Doppler_factor
      # Interpolate to original wavelength grid
      interpolationfunction = interp1d(
                new_wave, AbsorptionLine, kind="cubic", fill_value="extrapolate")
      interpolatedModel = interpolationfunction(wave)
      this_c, _ = pearsonr(flux, interpolatedModel)
      all_corr[loop] = this_c
      #print(v_rad, this_c)
      #plt.plot(wave, flux)
      #plt.gca().get_xaxis().get_major_formatter().set_useOffset(False)
      #plt.plot(wave, interpolatedModel, color="orange", marker="*")
      #plt.show()
   # Use the maximum maximum correlation as the initial guess.  
   v_rad_guess = v_rad_grid[np.argmax(all_corr)]
   print(v_rad_guess)
   #plt.plot(v_rad_grid, all_corr, marker='*')
   #plt.axvline(v_rad_guess, color='red')
   #plt.show()

   # Now let's look at step two. Make a model that more or less fits the strongest line, 
   # then calculate residuals and search again. 
   # These parameters produce an OK fit. 
   lambda0 = [3302.369, 3302.978]
   f = [9.21e-03, 4.60e-03]
   gamma = [6.280e7, 6.280e7]
   b = [2.0]
   N = [1e14]
   v_rad = [22.2]
   v_resolution = 4.00

   AbsorptionLine = voigt_absorption_line(
    wave,
    lambda0=lambda0,
    b=b,
    N=N,
    f=f,
    gamma=gamma,
    v_rad=v_rad,
    v_resolution=v_resolution,
    )
   # Calculate residuals
   ResidualFlux = flux-AbsorptionLine + 1

   plt.plot(wave, flux)
   plt.gca().get_xaxis().get_major_formatter().set_useOffset(False)
   plt.plot(wave, AbsorptionLine, color="orange", marker="*")
   plt.plot(wave,ResidualFlux, color="green")
   plt.show()
   
   
