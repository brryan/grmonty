import numpy as np
import argparse
import scipy.interpolate as itp
from math import log10
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description='Takes output points and values and evaluates and estimates spectra.')

fpoints="points.npy"
fvalues="values.npy"
points=np.load(fpoints)
values=np.load(fvalues)

parser.add_argument('mdotp',type=float,help='Accretion rate (in Medd).')
parser.add_argument('mp',type=float,help='Mass (in solar masses)')
parser.add_argument('thetap',type=float,help='Obersver angle (in degrees).')

args = parser.parse_args()
logmp=log10(args.mp)
logmdotp=log10(args.mdotp)
thetap=args.thetap

thetai=0.3
thetaf=100
ntheta=6
theta=np.linspace(thetai,thetaf,num=ntheta)
print theta

nui=1.2355897e8
nuf=4.98969853e29
nnu=200
nu=np.logspace(log10(nui),log10(nuf),nnu)
lognu=map(log10,nu)

lognulnu=[itp.interpn(points,values,(logmp,logmdotp,thetap,i),method="linear",bounds_error=True)[0] for i in lognu]
nulnu=[pow(10,i) for i in lognulnu]
ax = plt.subplot(1,1,1)
ax.set_title(r'$\dot{M}=10^{%1.3f}\dot{M}_{Edd}$,   $M=10^{%1.3f}M_{Sun}$,  $\theta=%1.2f\degree$' %(logmdotp, logmp, thetap))
ax.step(nu,nulnu, where='mid', color='k', linewidth=1)
ax.loglog()
ax.set_xlim([1.e8, 1.e22])
ax.set_ylim([1.e25, 1.e40])
ax.set_xlabel(r'$\nu$'+' (Hz)')
ax.set_ylabel(r'$\nu L\nu$'+' (erg s'+r'$^{-1}$'+')')
plt.show()