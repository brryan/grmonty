import numpy as np
import argparse
import scipy.interpolate as itp
from math import log10
import matplotlib.pyplot as plt
import subprocess

parser = argparse.ArgumentParser(description='Takes output points and values, estimates spectra, and compares it to grmonty output.')

fpoints="points.npy"
fvalues="values.npy"
points=np.load(fpoints)
values=np.load(fvalues)

parser.add_argument('dumpnam',type=str,help='Dump file location.')
parser.add_argument('mdotp',type=float,help='Accretion rate (in Medd).')
parser.add_argument('mp',type=float,help='Mass (in solar masses)')
parser.add_argument('thetabin',type=int,help='Theta bin.')

args = parser.parse_args()
dumpnam=args.dumpnam
logmp=log10(args.mp)
logmdotp=log10(args.mdotp)
thetabin=args.thetabin

nph=1e5
nui=1.2355897e8
nuf=4.98969853e29
nnu=200
nu=np.logspace(log10(nui),log10(nuf),nnu)
lognu=map(log10,nu)

thetai=0.3
thetaf=100
ntheta=6
theta=np.linspace(thetai,thetaf,num=ntheta)

lognulnu=[itp.interpn(points,values,(logmp,logmdotp,theta[thetabin],i),method="linear",bounds_error=True)[0] for i in lognu]
nulnu=[pow(10,i) for i in lognulnu]

#Run grmonty and import data in the exact same way the other plot program does.
subprocess.call(["./grmonty", repr(nph), dumpnam, repr(args.mdotp), repr(args.mp), repr(1)])

ME  = 9.1093897e-28
CL  = 2.99792458e10
HPL = 6.6260755e-27
LSUN = 3.827e33
N_THBINS = 6
BIN=thetabin

fnam = 'grmonty.spec'
data = np.loadtxt(fnam)
originalnu = 10.**(data[:,0])*ME*CL**2/HPL
DNULNU = (len(data[0])-1)/N_THBINS # File format offset between theta bins
originalnulnu = data[:,1 + BIN*DNULNU]*LSUN # erg s^-1

diff=[abs(nulnu[i]-originalnulnu[i])/originalnulnu[i] for i in range(0,nnu)]

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
ax.set_xlabel(r'$\nu$'+' (Hz)')
ax.set_ylabel(r'$\nu L\nu$'+' (erg s'+r'$^{-1}$'+')')
ax=plt.subplot(212)
ax.step(nu,diff, where='mid', color='k', linewidth=1)
ax.loglog()
ax.set_xlim([1.e8, 1.e22])
ax.set_xlabel(r'$\nu$'+' (Hz)')
#ax.set_ylabel(r'$(\nu L\nu_{Estimated}-\nu L\nu_{Calculated})/\nu L\nu_{Calculated}$')
ax.set_ylabel('Relative error')
plt.show()