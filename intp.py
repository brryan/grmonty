import numpy as np
import argparse
import scipy.interpolate as itp

parser = argparse.ArgumentParser(description='Takes output points and values and evaluates an interpolation at the requested point.')

fpoints="points.npy"
fvalues="values.npy"
points=np.load(fpoints,allow_pickle=False)
values=np.load(fvalues,allow_pickle=False)

parser.add_argument('mp',type=float,help='Mass (in solar masses)')
parser.add_argument('mdotp',type=float,help='Accretion rate (in Medd).')
parser.add_argument('thetap',type=float,help='Obersver angle (in degrees).')
parser.add_argument('nup',type=float,help='Frequency.')
args = parser.parse_args()
mp=args.mp
mdotp=args.mdotp
thetap=args.thetap
nup=args.nup

print itp.interpn(points,values,(mp,mdotp,thetap,nup),method="linear",bounds_error=True)