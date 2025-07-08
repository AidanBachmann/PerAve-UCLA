# Description: Main function for 1D FEL simulation. Adapted from Matlab version.

# ---------- Imports ----------

import Params as params # Parameters for code
import Functions as funcs # Function definitions
import numpy as np
from numba import jit
import time

# ---------- Main Script ----------

# Note: res_phase and Kz are numpy arrays (which are objects) so, by default, Python passes by reference. As a result,
# we don't need any return statements on functions such as computeUndulatorField.

## Flags
plotFilter = False # Set to true to plot filter function
npasses = int(100) # Number of passes through oscillator

## Arrays to store tapering of params

res_phase = np.zeros([params.Nsnap]) # Resonant phase?
Kz = np.zeros([params.Nsnap]) # Number of buckets for spatial discretization, I think

funcs.computeUndulatorField(res_phase,Kz) # Compute undulator field

## FEL Parameters

rho1D = (1/params.gamma0)*pow( (1/8)*(params.I/params.IA)*(pow(params.K,2)/(pow(params.sigmax*params.ku,2))),(1/3))
Lgain = params.lambdau/(4*np.sqrt(3)*np.pi*rho1D)
Lsat =   params.lambdau/rho1D
Psat = 1.6*rho1D*params.Ee*params.I
if params.tapering != int(0):
    bucket_params = funcs.bucket_parameters()
    # *** Insert code for other parameters here ***
Lgain3D = funcs.correction3D(Lgain) # Compute 3D correction to Lgain

# *** 3D correction to Lgain goes here ***

## Run the main integration routine

cavitydetuning = -16 # In units of zsep
transmission = 0.66 # Power transmission through one cavity pass, losses = 1 - transmission                                
sigma_omega = 0.003*params.nslices*params.zsep # Filter fractional bandwidth

filter3 = funcs.filter(sigma_omega,plotFilter) # Compute filter

funcs.oscLoop(npasses,Kz,res_phase) # Oscillator loop