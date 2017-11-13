import sys
import numpy as np
import argparse
import matplotlib.pyplot as plt
import math
import subprocess

dumpnam="../grmhd_data/dump040"
fnam="grmonty.spec"

# CGS constants
LSUN = 3.827e33

# Take input arguments
parser = argparse.ArgumentParser(description='Runs grmonty on a grid of the input variables and outputs a set of points and values.')

## Define input variables
parser.add_argument('nph',type=float,help='Number of photons to use.')
# Accretion rates
parser.add_argument('mdoti',type=float,help='Lower accretion limit (in Medd).')
parser.add_argument('mdotf',type=float,help='Upper accretion limit (in Medd).')
parser.add_argument('nmdot',type=float,help='Amount of points to use (including initial and final).')
# BH masses
parser.add_argument('mi',type=float,help='Lower mass limit (in Msun).')
parser.add_argument('mf',type=float,help='Upper mass limit (in Msun)n.')
parser.add_argument('nm',type=float,help='Amount of points to use (including initial and final).')

## Take input variables
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

#nui=1e8
#nuf=1e30
nnu=200

#Create vectors for all variables
mdot=np.logspace(math.log(mdoti,10),math.log(mdotf,10),nmdot,endpoint=True)
m=np.logspace(math.log(mi,10),math.log(mf,10),nm,endpoint=True)
theta=np.linspace(thetai,thetaf,num=ntheta,endpoint=True)
nu=np.logspace(math.log(1.2355897e8,10),math.log(4.98969853e29,10),nnu,endpoint=True)

# Declare arrays
specdata=np.empty((len(m),len(mdot),ntheta,nnu),dtype=float)
points=[]
values=[]

# Logarithm of the vectors for interpolation
logm=[math.log10(i) for i in m]
logmdot=[math.log10(i) for i in mdot]
lognu=[math.log10(i) for i in nu]

#Call grmonty iteratively.
for i in range(0,len(m)):
    for j in range(0,len(mdot)):
        subprocess.call(["./grmonty", repr(nph), dumpnam, repr(mdot[j]), repr(m[i]), repr(1)])
        data=np.loadtxt(fnam)
        DNULNU = (len(data[0])-1)/ntheta
        for k in range(0,len(theta)):
            specdata[i][j][k][:]=data[:,1 + k*DNULNU]*LSUN
            for l in range(0,nnu):
                points.append((logm[i],logmdot[j],theta[k],lognu[l]))
                if specdata[i][j][k][l]>0:
                    values.append(math.log10(specdata[i][j][k][l]))

# Export results
np.save("points",points,allow_pickle=False)
np.save("values",values,allow_pickle=False)