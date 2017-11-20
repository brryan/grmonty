import numpy as np
from math import log10
from math import sqrt
import matplotlib.pyplot as plt
import subprocess
import argparse

# Physical constants
ME  = 9.1093897e-28
CL  = 2.99792458e10
HPL = 6.6260755e-27
LSUN = 3.827e33

# grmonty constants
ntheta = 6
nnu=200
fnam = 'grmonty.spec'

# Input
parser = argparse.ArgumentParser(description='Run grmonty with different numbers of superphotons and check N^(-1/2) convergence.')

parser.add_argument('dumpnam',type=str,help='Dump file location.')
parser.add_argument('ni',type=float,help='Initial number of superphotons.')
parser.add_argument('nf',type=float,help='Final number of superphotons.')
parser.add_argument('nn',type=int,help='Amount of different superphoton numbers to use.')
parser.add_argument('mdotp',type=float,help='Accretion rate (in Medd).')
parser.add_argument('mp',type=float,help='Mass (in solar masses)')
parser.add_argument('thetabin',type=int,help='Theta bin to compare.')

args = parser.parse_args()
dumpnam=args.dumpnam
ni=args.ni
nf=args.nf
nn=args.nn
logmp=log10(args.mp)
logmdotp=log10(args.mdotp)
thetabin=args.thetabin
nph=args.nph

# Super photon numbers to compare
ns=np.logspace(log10(ni),log10(nf),nn)

# Declare array that contains all spectrum data: nulnu[run number (corresponding to number of superphotons used),frequency]
nulnu=np.empty((nn,nnu),dtype=float)

# L1 norm vector
l=np.zeros(shape=(nn-1,))

# Call grmonty for every number of superphotons
for i in range(0,nn):
    subprocess.call(["./grmonty", repr(ns[i]), dumpnam, repr(mdot), repr(m), repr(1)])
    data = np.loadtxt(fnam)
    nu = 10.**(data[:,0])*ME*CL**2/HPL
    DNULNU = (len(data[0])-1)/ntheta # File format offset between theta bins
    # Store in corresponding vector
    nulnu[i][:] = data[:,1 + thetabin*DNULNU]*LSUN

# Compute L1 norm
for i in range(1,nn):
    l[i-1]=np.sum([abs(nulnu[i][j]-nulnu[i-1][j]) for j in range(0,nnu)])

# Plot N^(-1/2)
c=l[0]/(ns[1]**(-0.5))
def f(t):
    return c*(t**(-0.5))

y=map(f,ns[1:nn])

# Plot number of superphotons and L1 norm
plt.title('grmonty convergence')
plt.plot(ns[1:nn],l,'ko')
plt.plot(ns[1:nn],y,'--k',label=r'$N^{-1/2}$')
plt.legend(loc="upper right")
plt.loglog()
ax=plt.gca()
ax.set_xlabel(r'$N$')
ax.set_ylabel(r'$L_{1}$')
plt.show()