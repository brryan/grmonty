import numpy as np
import argparse
from math import log10
import os
import subprocess

import matplotlib.pyplot as plt

fnam="spectrum.dat"

# CGS constants
LSUN = 3.827e33

# Input
parser = argparse.ArgumentParser(description='Runs grmonty on a grid of the input variables and outputs a set of support points for interpolation.')

parser.add_argument('nph',type=float,help='Number of photons to use.')
parser.add_argument('dumpnam',type=str,help='Dump file location.')

parser.add_argument('mdoti',type=float,help='Lower accretion limit (in Medd).')
parser.add_argument('mdotf',type=float,help='Upper accretion limit (in Medd).')
parser.add_argument('nmdot',type=int,help='Amount of points to use (including initial and final).')

parser.add_argument('mi',type=float,help='Lower mass limit (in Msun).')
parser.add_argument('mf',type=float,help='Upper mass limit (in Msun).')
parser.add_argument('nm',type=int,help='Amount of points to use (including initial and final).')

args = parser.parse_args()
nph=args.nph
dumpnam=args.dumpnam
mdoti=args.mdoti
mdotf=args.mdotf
nmdot=args.nmdot
mi=args.mi
mf=args.mf
nm=args.nm

# Observer angles
thetai=0.3
thetaf=100
ntheta=6
NVAR=8
# Frequencies
nui=1e8
nuf=1e24
nnu=200

#Create vectors for all variables,
mdot=np.logspace(log10(mdoti),log10(mdotf),nmdot)
m=np.logspace(log10(mi),log10(mf),nm)
theta=np.linspace(thetai,thetaf,ntheta)
nu=np.logspace(log10(nui),log10(nuf),nnu)
# Logarithm of the vectors for interpolation
logm=map(log10,m)
logmdot=map(log10,mdot)
lognu=map(log10,nu)

# Declare arrays containing spectrum information
values=np.empty((nm,nmdot,ntheta,nnu),dtype=float)
points=[logm,logmdot,theta,lognu]

#Call grmonty iteratively and store results.
FNULL = open(os.devnull, 'w')
for i in range(0,nm):
    for j in range(0,nmdot):
        print("Currently running grmonty with N="+repr(nph)+", Mdot="+repr(mdot[j])+", M="+repr(m[i])+".")
        subprocess.call(["./grmonty", repr(nph), dumpnam, repr(mdot[j]), repr(m[i]), repr(1)],stdout=FNULL, stderr=subprocess.STDOUT)
        data=np.loadtxt(fnam)
        nbin = (len(data[0])-1)/NVAR
        for k in range(0,ntheta):
            values[i][j][k][:]=data[:,1 + k*NVAR]*LSUN
            for l in range(0,nnu):
                if values[i][j][k][l]>0:
                    values[i][j][k][l]=log10(values[i][j][k][l])

# Export results
np.save("points",points)
np.save("values",values)

#plt.step(nu,[pow(10,values[nm-1][nmdot-1][1][i]) for i in range(0,nnu)],'k')
#plt.loglog()
#ax=plt.gca()
#ax.set_xlim([1.e8, 1.e22])
#ax.set_ylim([1.e25, 1.e40])
#plt.show()
