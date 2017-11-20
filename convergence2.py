import argparse
from math import log10
from math import sqrt
import matplotlib.pyplot as plt
import numpy as np
import subprocess

# Physical constants
ME  = 9.1093897e-28
CL  = 2.99792458e10
HPL = 6.6260755e-27
LSUN = 3.827e33

# Dump and output file names
fnam = 'grmonty.spec'

# Input
parser = argparse.ArgumentParser(description='Compare two grmonty runs for each number of superphotons.')

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
mp=args.mp
mdotp=args.mdotp
thetabin=args.thetabin

# Theta bins and number of frequencies
N_THBINS = 6
nnu=200

# Super photon numbers to compare
ns=np.logspace(log10(ni),log10(nf),nn)

# Declare array that contains all spectrum data: nulnu[run number (corresponding to number of superphotons used),frequency]
nulnu=np.empty((nn,2,nnu),dtype=float)
l=np.zeros(shape=(nn,))

#Call grmonty for every number of superphotons
for i in range(0,nn):
    subprocess.call(["./grmonty", repr(ns[i]), dumpnam, repr(mdotp), repr(mp), repr(1)])
    data = np.loadtxt(fnam)
    nu = 10.**(data[:,0])*ME*CL**2/HPL
    DNULNU = (len(data[0])-1)/N_THBINS # File format offset between theta bins
    # Store in corresponding vector
    nulnu[i][0][:] = data[:,1 + thetabin*DNULNU]*LSUN
    subprocess.call(["./grmonty", repr(ns[i]), dumpnam, repr(mdotp), repr(mp), repr(1)])
    data = np.loadtxt(fnam)
    nu = 10.**(data[:,0])*ME*CL**2/HPL
    DNULNU = (len(data[0])-1)/N_THBINS # File format offset between theta bins
    nulnu[i][1][:] = data[:,1 + thetabin*DNULNU]*LSUN

# Compute L1 norm
for i in range(0,nn):
    l[i]=np.sum([abs(nulnu[i][0][j]-nulnu[i][1][j]) for j in range(0,nnu)])

# Plot N^(-1/2)
c=l[0]/(ns[0]**(-0.5))
def f(t):
    return c*(t**(-0.5))

x=np.logspace(log10(ni),log10(nf),100)
y=map(f,x)
print l
# Plot number of superphotons and L1 norm
plt.plot(ns,l,'ko')
plt.plot(x,y,'k--')
plt.loglog()
#plt.plot(x,y2,'--k')
ax=plt.gca()
ax.set_xlabel(r'$N$')
ax.set_ylabel(r'$\sum\Delta\nu L\nu$')
plt.show()