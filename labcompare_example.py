from edibles.edibles import DATADIR
from edibles.edibles import PYTHONDIR
from edibles.edibles.utils.edibles_spectrum import EdiblesSpectrum
from edibles.edibles.utils.edibles_oracle import EdiblesOracle
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from PyAstronomy import pyasl

# Example: Look at Farid's spectrum of pentacene. 
labfilename = PYTHONDIR + '/edibles/data/Labdata/CRDS/PENTACENE.DAT'
# The file contains wavenumbers (in vacuum) and intensities. Read those in as pandas dataframe. 
labspec = pd.read_csv(labfilename, delim_whitespace=True)
# Add a column to contain the (air) wavelength in AA. 
labspec['wavelength'] = pyasl.vactoair2(1e8/labspec['wno'],mode='ciddor')
normint = labspec.int - np.median(labspec.int)
labspec['norm'] = normint / np.max(normint)

# Pentacene has strongest bands at about 5338 and 5361 AA. Let's search what EDIBLES spectra
# We have in that wavelength range and compare. 
plotrange = [5330,5380]
pythia = EdiblesOracle()
List = pythia.GetObsListByWavelength(5338)
for filename in List: 
	sp = EdiblesSpectrum(filename)
	wave = sp.wave
	flux = np.clip(sp.flux, 0, None) 
	bool_keep = (wave > plotrange[0]) & (wave < plotrange[1])
	plotwave = wave[bool_keep]
	plotflux = flux[bool_keep]
	normflux = plotflux / np.median(plotflux)
	plt.figure()
	plt.xlim(5330,5380)
	ylim = [np.min(normflux), np.max(normflux)]
	plt.ylim(ylim)
	plt.plot(plotwave, normflux)
	# Rescale lab spectrum to plot range
	dynrange = ylim[1]-ylim[0]
	plt.plot(labspec.wavelength, labspec.norm * dynrange/5 + 1)
	plt.show()

#plt.plot(labspec.wavelength,labspec.int)
#plt.show()