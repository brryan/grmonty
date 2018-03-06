from hdf5_to_dict import *
from math import pi,log10
import numpy as np
from units import *
import matplotlib.pyplot as plt

# Setup
fname = 'dump_00002000.h5'
hdr = load_hdr(fname)
geom = load_geom(hdr, recalc = False)
gdet = geom['gdet']
dump = load_dump(fname,geom)
cgs = get_cgs()

# Parameters from the grmonty run
TPTE = 3
RHO_unit = 8.495121e-17
L_unit = 6.646685e+11

# Frequencies
nu =  np.logspace(log10(1.23558978e+08),log10(4.98969853e+29),200)

# Units
Ne_unit = RHO_unit / cgs['MP']
Thetae_unit = ( ( hdr['gam'] - 1 ) * cgs['MP'] / cgs['ME'] ) / ( 1 + TPTE )
print 'Thetae_unit = ',Thetae_unit,' Ne_unit = ',Ne_unit

# Number density and Te
Ne = dump['RHO'] * Ne_unit
Thetae = dump['UU'] / dump['RHO'] * Thetae_unit
Te = Thetae * cgs['ME'] * cgs['CL'] * cgs['CL'] / cgs['KBOL']

# Calculation done in jnu_bremss in jnu_mixed.c
efac = np.array([np.exp(- cgs['HPL'] * i / (cgs['KBOL'] * Te)) for i in nu])
rel = 1.0 + 4.4e-10 * Te
gff = 1.2
jv0 = 32 * pi * np.power(cgs['QE'],6) / (12.0 * pi * cgs['ME'] * np.power(cgs['CL'],3))
jv0 *= np.power(2*pi/(3*cgs['KBOL']*cgs['ME']),0.5)
jv0 *= np.power(Te,-0.5) * Ne * Ne
jv = np.array([efac[i,:,:,:] * rel * gff * jv0 for i in range(0,len(nu))])

# Integration
integral = np.empty(len(nu))
for i in range(0,len(nu)):
	integral[i] = nu[i] * 4 * pi * hdr['dx1'] * hdr['dx2'] * hdr['dx3'] * ( L_unit ** 3 ) * (jv[i,:,:,:] * gdet[:,:,None]).sum()

np.save('BremssEstimateTPTE'+str(TPTE)+'.dat',integral)
#ax = plt.gca()
#ax.step(nu, integral, where='mid', color='k', linewidth=1)
#ax.set_xscale('log')
#ax.set_yscale('log')
#ax.set_xlim([1.e8, 1.e24])
#ax.set_ylim([1.e25, 1.e40])
#ax.set_xlabel('nu (Hz)')
#ax.set_ylabel('nuLnu (erg s^-1)')
#plt.show()
