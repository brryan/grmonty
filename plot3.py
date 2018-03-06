import numpy as np
import matplotlib.pyplot as plt

ME   = 9.1093897e-28
CL   = 2.99792458e10
HPL  = 6.6260755e-27
LSUN = 3.827e33

fnam1 = 'spectrum.dat'
fnam2 = 'BremssEstimateTPTE10.dat'

data1 = np.loadtxt(fnam1)
data2 = np.loadtxt(fnam2)
nu = 10.**(data1[:,0])*ME*CL**2/HPL

NVAR = 8
nbin = (len(data1[0])-1)/NVAR
error=np.empty(len(nu))
nuLnu1 = np.array(data1[:,NVAR+1]*LSUN)
nuLnu2 = np.array(data2)

for i in range(0,len(nu)):
	if nuLnu1[i] != 0:
		error[i]=abs(nuLnu2[i]-nuLnu1[i])/nuLnu1[i]
	else:
		error[i]=1e-5

ax = plt.subplot(3,1,1)
ax.set_title(r'$N=10^{6}$, $M=4.5\times10^{6}M_{Sun}$, $\dot{M}=10^{-8}\dot{M}_{Edd}$')
ax.step(nu, nuLnu1, where='mid',color='k')
ax.set_xscale('log'); ax.set_yscale('log')
nuLnu_max = nuLnu1.max()
ax.set_ylim([1.e-10*nuLnu_max, 1.e1*nuLnu_max])
ax.set_xlim([1.e8, 1.e24])
ax.set_ylabel(r'$\nu L\nu$')
ax2 = plt.subplot(3,1,2)
ax2.step(nu,nuLnu2,where='mid',color='k')
ax2.loglog()
ax2.set_ylabel(r'$\nu L\nu$')
ax2.set_ylim([1e-10*nuLnu_max,10*nuLnu_max])
ax2.set_xlim([1e8,1e24])
ax3=plt.subplot(3,1,3)
ax3.set_ylabel(r'$\epsilon$')
ax3.set_xlabel(r'$\nu$')
ax3.step(nu,error,where='mid',color='k')
ax3.set_ylim([1e-5,10*error.max()])
ax3.set_xscale('log')
ax3.set_yscale('log')
ax3.set_xlim([1e8,1e24])
plt.show()
plt.savefig('CompareBremssgGrmontyEstimateTPTE10.png')
