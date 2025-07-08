# Description: Define parameters for 1D FEL Code.

# ---------- Imports ----------

import numpy as np
import math

# ---------- Physical Constants ----------

c = 2.99792458e8 # Speed of light
e0 = 1.60217657e-19 # Elementary charge
me = 9.10938291e-31 # Electron mass
eps0 = 8.85418782e-12 # Permittivity of free space
mu0 = 1.256637e-6 # Permeability of free space
IA = 17045 # Alfven current
Z0 = 376.73 # Impedance of free space

# ---------- Parameters ----------

sigma_t = 0.9e-12 # Pulse width
use3Dcorrection  = int(1) # Use 3D correction boolean
beamdistribution = int(1) # Set beam distribution, 2-uniform 1-gaussian
laserdistribution = int(1) # Set beam distribution, 2-uniform 1-gaussian
recirculate = int(0) # ***** Not sure

## Undulator parameters
lambdau = 3.14e-2 # Undulator period (m)
K = 1.934 # RMS undulator parameter, e0*Bfield*/me/c/ku
ku = (2*np.pi)/lambdau # Undulator wavenumber
lwig = 9.42 # Undulator length (m), paper defines number of undulator periods and lambdau, so lwig is computed from these vals
# Tapering options
tapering = int(0) # Set tapering (-1 acceleration, 0 no tapering, 1 decelation)   
tapering_strength = int(2)   # 0 --> max of slices at time 0, 1 --> max of slices, 2 --> avg of slices 
z0 = 0
psir = np.pi/6 # Phase

## Simulation control options
phasespacemovie = int(0) # Create phase space movie
itdp = int(0) # 1 for time-dependent simulation, 0 for time independent
prebunching = int(1) # Prebunching boolean, set to 1 to start from a pre-bunched beam
changeresphase = int(0) # Update resonant phase
saveoutput = int(0) # Save output bool
# Set simulation length and # of snapshots
delz = 1 # ***** Not sure
und_periods = round(lwig/lambdau) # Number of undulator periods to simulate
Nsnap = round(lwig/lambdau/delz) # Number of snapshots to take over the length of the undulator
zsep = 5 # ***** Not sure
Nslip = round(Nsnap/zsep) # ***** Not sure
shotnoise = 1 # ***** Not sure
lambda0 = 13.5e-9 # Seed wavelength (m)
k = (2*np.pi)/lambda0 # Wavenumber in free space
nslices = 4*Nslip + 4*round((sigma_t/(zsep*lambda0))*c) # ***** Not sure

if itdp == 0:
    print('\nRunning time independent simulation.\n')
    nslices = 1
elif itdp == 1:
    print('\nRunning time-dependent simulation.\n')
    Nsnap = math.floor(und_periods/delz) # Note that the result must be an integer
stepsize = lambdau*delz

## Electron beam parameters
gamma0 = np.sqrt((k/(2*ku))*(1+pow(K,2))) # Lorentz factor
Np = int(512) # Number of macroparticles (500-1000 well)
Ee = (gamma0*me*pow(c,2))/e0 # Total e-beam energy, gamma*m*c^2 (eV)
energyspread = 1*20e-15/sigma_t # Absolute energy spread (MeV)
deltagammarel = energyspread/(gamma0*0.511) # Relative energy spread, dgamma/gamma
deltagamma = gamma0*deltagammarel
bunch = 0.7 # Initial bunching factor
bunchphase = -psir-np.pi/2 # Initial bunching phase
buncherAmp = 5 # ***** Not sure

betax = 5.1 # Average beta function, 4.7544975 (something weird going on here)
emitx = 1e-6 # normalized emittance in m*rad
charge = 2255.965e-12 # Computed from peak current and sigma_t
if (beamdistribution == 1): # Set current based on beam distribution
    I = charge/(np.sqrt(2*np.pi)*sigma_t) # Beam current 
else:
    I = charge/(2*sigma_t) # Beam current

sigmax = np.sqrt(betax*emitx/gamma0) # Beam radius
A_e = 2*np.pi*pow(sigmax,2) # Beam cross section 
Simulation_temporal_window = nslices*zsep*(lambda0/c)

## Radiation parameters
P0 = 30e9 # Peak input power (W)
A_mode = A_e # 1D code, area is same for e_beam and radiation
waist = np.sqrt((2*A_mode)/np.pi) # Waist size
zr = (np.pi*pow(waist,2))/lambda0 # Rayleigh length of seed (m)
E0 = np.sqrt((2*P0)/(2*c*eps0*A_mode)) # Assume circular polarization
slippage = (nslices/2)*( (lambda0*zsep)/c ) # ***** Not sure
sigma_l = 2400e-15 # ***** Not sure

## Simplifying constants
chi2 = e0/(me*pow(c,2))
chi1 = ((mu0*c)/2)*(I/A_e)


## Initializing simulation arrays based on defined parameters
# All arrays defined below were previously computed in each pass of PerAve core. It is much more efficient to compute them once and resuse the arrays.

n_electron = (I*lambda0*zsep)/(e0*c) # Number of electrons?
p1 = np.zeros([Np,1])

tslice = np.linspace(1,nslices,nslices,dtype='int')*( (lambda0*zsep)/c ) # Time slices

if (beamdistribution == 1): # Compute beam distribution
    profile_b = np.exp(-pow(tslice-tslice[-1]*np.ones(nslices)/2,2)/pow(2*sigma_t,2))
else:
    profile_b = np.zeros([nslices])
    profile_b[abs(tslice-tslice[-1]*np.ones([nslices])/2)<sigma_t] = 1

if (laserdistribution == 1): # Compute laser distribution
    profile_l = np.exp(-pow(tslice-slippage*np.ones([nslices]),2)/(2*pow(sigma_l,2)))
else:
    profile_l = np.zeros([nslices])
    profile_l[abs(tslice-slippage*np.ones([nslices]))<sigma_l] = 1