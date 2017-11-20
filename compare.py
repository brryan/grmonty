import numpy as np
import argparse
import scipy.interpolate as itp
from math import log10
import matplotlib.pyplot as plt
import subprocess

# Physical constants
ME  = 9.1093897e-28
CL  = 2.99792458e10
HPL = 6.6260755e-27
LSUN = 3.827e33

# Files where support points are stored
points=np.load("points.npy")
values=np.load("values.npy")
# grmonty output file
fnam = 'grmonty.spec'

# Input
parser = argparse.ArgumentParser(description='Takes output points and values, estimates spectra, and compares it to grmonty output.')

parser.add_argument('dumpnam',type=str,help='Dump file location.')
parser.add_argument('nph',type=int,help='Number of photons to use in grmonty run.')
parser.add_argument('mdotp',type=float,help='Accretion rate (in Medd).')
parser.add_argument('mp',type=float,help='Mass (in solar masses)')
parser.add_argument('thetabin',type=int,help='Theta bin to compare.')

args = parser.parse_args()
dumpnam=args.dumpnam
mp=args.mp
mdotp=args.mdotp
thetabin=args.thetabin
nph=args.nph

logmp=log10(mp)
logmdotp=log10(mdotp)

# Frequencies vector
nui=1.2355897e8
nuf=4.98969853e29
nnu=200
nu=np.logspace(log10(nui),log10(nuf),nnu)
lognu=map(log10,nu)

# Viewing angles vector
thetai=0.3
thetaf=100
ntheta=6
theta=np.linspace(thetai,thetaf,ntheta)

# Interpolated values for nuLnu
lognulnu=[itp.interpn(points,values,(logmp,logmdotp,theta[thetabin],i),method="linear",bounds_error=True)[0] for i in lognu]
nulnu=[pow(10,i) for i in lognulnu]

#Run grmonty
subprocess.call(["./grmonty", repr(nph), dumpnam, repr(mdotp), repr(mp), repr(1)])

# nuLnu for the grmonty run
data = np.loadtxt(fnam)
originalnu = 10.**(data[:,0])*ME*CL**2/HPL
DNULNU = (len(data[0])-1)/ntheta # File format offset between theta bins
originalnulnu = data[:,1 + thetabin*DNULNU]*LSUN # erg s^-1

# Relative error between interpolated and grmonty outputs
err=[abs(nulnu[i]/originalnulnu[i]-1) for i in range(0,nnu)]
# Discard unreasonably high error
for i in range(0,nnu):
    if err[i]>1e2:
        err[i]=0

# Average relative error
averr=[np.average(err) for i in range(0,nnu)]

plt.figure(1)
ax=plt.subplot(211)
ax.set_title(r'$N_{ph}=10^{%1.2f}$,  $\dot{M}=10^{%1.2f}\dot{M}_{Edd}$,  $M=10^{%1.2f}M_{Sun}$,  $\theta=%1.2f\degree$' %(log10(nph),logmdotp, logmp, theta[thetabin]))
ax.step(nu,nulnu, where='mid', color='k', linewidth=1,label='Interpolated')
ax.legend(loc="upper right")
ax.step(originalnu,originalnulnu, where='mid', color='g', linewidth=1,label='grmonty')
ax.legend(loc="upper right")
ax.loglog()
ax.set_xlim([1.e8, 1.e22])
ax.set_ylim([1.e25, 1.e40])
ax.set_xlabel(r'$\nu$ (Hz)')
ax.set_ylabel(r'$\nu L\nu$ (erg s$^{-1}$)')
ax=plt.subplot(212)
ax.step(nu,err, where='mid', color='k', linewidth=1)
ax.plot(nu,averr,'--', color='k', linewidth=1,label='Average')
ax.legend(loc="upper right")
ax.loglog()
ax.set_xlim([1.e8, 1.e22])
ax.set_xlabel(r'$\nu$ (Hz)')
#ax.set_ylabel(r'$(\nu L\nu_{Estimated}-\nu L\nu_{Calculated})/\nu L\nu_{Calculated}$')
ax.set_ylabel('Relative error')
plt.show()