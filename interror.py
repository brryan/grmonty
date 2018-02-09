import argparse
from math import log10
import matplotlib.pyplot as plt
import numpy as np
import os
import scipy.interpolate as itp
import subprocess

fnam="spectrum.dat"
FNULL = open(os.devnull, 'w')

# CGS constants
LSUN = 3.827e33
ME   = 9.1093897e-28
CL   = 2.99792458e10
HPL  = 6.6260755e-27

# Input
parser = argparse.ArgumentParser(description='Runs grmonty interpolation using a grid of decreasing size.')

parser.add_argument('nph',type=float,help='Number of photons to use.')
parser.add_argument('dumpnam',type=str,help='Dump file location.')

parser.add_argument('mdoti',type=float,help='Lower accretion limit (in Medd).')
parser.add_argument('mdotf',type=float,help='Upper accretion limit (in Medd).')

parser.add_argument('mi',type=float,help='Lower mass limit (in Msun).')
parser.add_argument('mf',type=float,help='Upper mass limit (in Msun).')
parser.add_argument('nit',type=int,help='Number of iterations.')

parser.add_argument('mdotp',type=float,help='Accretion rate (in Medd).')
parser.add_argument('mp',type=float,help='Mass (in solar masses)')
parser.add_argument('thetabin',type=int,help='Theta bin.')

args = parser.parse_args()
nph=args.nph
dumpnam=args.dumpnam
mdoti=args.mdoti
mdotf=args.mdotf
mi=args.mi
mf=args.mf
nit=args.nit
mp=args.mp
logmp=log10(mp)
mdotp=args.mdotp
logmdotp=log10(mdotp)

thetabin=args.thetabin

# Observer angles
thetai=0.3
thetaf=100
ntheta=6
NVAR=8
theta=np.linspace(thetai,thetaf,ntheta)
thetap=theta[thetabin]

# Run grmonty on the point of interest twice and get the error
print("Currently running grmonty with N="+repr(nph)+", Mdot="+repr(mdotp)+", M="+repr(mp)+".")
subprocess.call(["./grmonty", repr(nph), dumpnam, repr(mdotp), repr(mp), repr(1)],stdout=FNULL, stderr=subprocess.STDOUT)
data=np.loadtxt(fnam)
# Frequencies
#nu = 10.**(data[:,0])*ME*CL**2/HPL
#nnu=len(nu)
nui=1e8
nuf=1e24
nnu=200
nu=np.logspace(log(nui),log(nuf),nnu)
lognu=map(log10,nu)

nbin = (len(data[0])-1)/NVAR
nulnu1=data[:,1 + thetabin*NVAR]*LSUN
print("Currently running grmonty with N="+repr(nph)+", Mdot="+repr(mdotp)+", M="+repr(mp)+".")
subprocess.call(["./grmonty", repr(nph), dumpnam, repr(mdotp), repr(mp), repr(1)],stdout=FNULL, stderr=subprocess.STDOUT)
data=np.loadtxt(fnam)
nbin = (len(data[0])-1)/NVAR
nulnu2=data[:,1 + thetabin*NVAR]*LSUN
baseerror=np.array([abs(nulnu1[i]-nulnu2[i]) for i in range(0,nnu)])

#Create vectors for all variables,
mdot=np.logspace(log10(mdoti),log10(mdotf),2*nit)
m=np.logspace(log10(mi),log10(mf),2*nit)

errors=np.empty((nit,85))
lognulnuint=np.empty((nit,nnu))
nulnuint=np.empty((nit,nnu))
averrors=np.empty(nit)

for i in range(0,nit):
    print("Current grid: Mdoti= "+repr(mdot[i])+" , Mdotf="+repr(mdot[-i-1])+" , Mi= "+repr(m[i])+", Mf="+repr(m[-i-1]))
    subprocess.call(["python", "iterate.py", repr(nph), dumpnam, repr(mdot[i]),repr(mdot[-i-1]),repr(2), repr(m[i]),repr(m[-i-1]),repr(2)])
    # Files containing grid and values
    points=np.load("points.npy")
    values=np.load("values.npy")
    lognulnuint[i,:]=np.array([itp.interpn(points,values,(logmp,logmdotp,thetap,j),method="linear",bounds_error=True)[0] for j in lognu])
    nulnuint[i,:]=np.array([pow(10,lognulnuint[i,j]) for j in range(0,nnu)])
    for j in range(18,nnu-97):
        if nulnu2[j]>1e27 and nulnuint[i,j]>1e27:
    	    errors[i,j-18]=abs(abs(nulnuint[i,j]-nulnu2[j])-abs(baseerror[j]))/nulnu2[j]
        else:
            errors[i,j-18]=-1
	print 'Error at nu=',j,": ",errors[i,j-18]
	if errors[i,j-18]>10:
		print 'HERE!'
    averrors[i]=np.mean(errors[i,:])
#averrors=np.nanmean(np.where(errors>=0,errors,np.nan),axis=1)
# Plot iteration number and average error
plt.title('Interpolation error')
plt.plot(range(0,nit),averrors,'k-')
ax=plt.gca()
ax.set_xlabel(r'$N_{Iteration}$')
ax.set_ylabel(r'$Average relative error$')
plt.show()

for i in range(0,nit):
	plt.figure(1)
	ax0=plt.subplot(211)
	ax0.set_title(r'$\epsilon=%1.3f$, Grid: $\dot{M}=[%1.2e,%1.2e], M=[%1.2e,%1.2e]$' %(averrors[i],mdot[i],mdot[-i-1],m[i],m[-i-1]))
	ax0.set_xscale('log')
	ax0.set_yscale('log')
	ax0.set_xlim([1.e8, 1.e24])
	ax0.set_ylim([1.e28, 1.e35])
	ax0.step(nu,nulnu2,where='mid', color='g', linewidth=1,label='grmonty')
	ax0.step(nu,nulnuint[i], where='mid', color='k', linewidth=1,label='Interpolated')
	ax0.legend(loc="upper right")
	ax0.set_xlabel(r'$\nu$')
	ax0.set_ylabel(r'$\nu L\nu$')

	ax1=plt.subplot(212)
#	ax1.set_title(r'$\epsilon$')
	ax1.set_xscale('log')
	ax1.set_yscale('log')
	ax1.set_xlim([1.e8, 1.e24])
	ax1.set_ylim([1e-5, 100])
	ax1.set_xlabel(r'$\nu$ (Hz)')
	ax1.set_ylabel(r'$\epsilon$')
	ax1.step(nu[18:nnu-97],errors[i,:],where='mid',color='k',linewidth=1,label=r'$\epsilon$')
	ax1.plot(nu,np.array([averrors[i] for j in range(0,nnu)]),'--', color='k', linewidth=1,label=r'$\bar{\epsilon}$')
	ax1.legend(loc="upper right")

	plt.savefig('Plot1'+str(i).zfill(3)+'.png')
	plt.show()
	plt.clf()
	plt.cla()
	plt.close()
