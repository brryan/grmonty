import numpy as np
import argparse
import matplotlib.pyplot as plt
from math import log10
import subprocess

dumpnam="../grmhd_data/dump040"
fnam="grmonty.spec"

# CGS constants
LSUN = 3.827e33

# Take input arguments
parser = argparse.ArgumentParser(description='Runs grmonty on a grid of the input variables and outputs a set of points and values.')

## Define input variables
parser.add_argument('nph',type=int,help='Number of photons to use.')
# Accretion rates
parser.add_argument('mdoti',type=float,help='Lower accretion limit (in Medd).')
parser.add_argument('mdotf',type=float,help='Upper accretion limit (in Medd).')
parser.add_argument('nmdot',type=int,help='Amount of points to use (including initial and final).')
# BH masses
parser.add_argument('mi',type=float,help='Lower mass limit (in Msun).')
parser.add_argument('mf',type=float,help='Upper mass limit (in Msun)n.')
parser.add_argument('nm',type=int,help='Amount of points to use (including initial and final).')

# Take input variables
args = parser.parse_args()
nph=args.nph
mdoti=args.mdoti
mdotf=args.mdotf
nmdot=args.nmdot
mi=args.mi
mf=args.mf
nm=args.nm

# Hardcoded in grmonty, will change in the near future.
thetai=0.3
thetaf=100
ntheta=6

nui=1.2355897e8
nuf=4.98969853e29
nnu=200

#Create vectors for all variables,
mdot=np.logspace(log10(mdoti),log10(mdotf),nmdot)
m=np.logspace(log10(mi),log10(mf),nm)
theta=np.linspace(thetai,thetaf,num=ntheta)
nu=np.logspace(log10(nui),log10(nuf),nnu)

# Logarithm of the vectors for interpolation
logm=[log10(i) for i in m]
logmdot=[log10(i) for i in mdot]
lognu=[log10(i) for i in nu]

# Declare arrays
values=np.empty((len(m),len(mdot),ntheta,nnu),dtype=float)
points=[logm,logmdot,theta,lognu]

#Call grmonty iteratively.
for i in range(0,len(m)):
    for j in range(0,len(mdot)):
        subprocess.call(["./grmonty", repr(nph), dumpnam, repr(mdot[j]), repr(m[i]), repr(1)])
        data=np.loadtxt(fnam)
        DNULNU = (len(data[0])-1)/ntheta
        for k in range(0,len(theta)):
            values[i][j][k][:]=data[:,1 + k*DNULNU]*LSUN
            for l in range(0,nnu):
                if values[i][j][k][l]>0:
                    values[i][j][k][l]=log10(values[i][j][k][l])

# Export results
np.save("points",points)
np.save("values",values)