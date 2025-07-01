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

res_phase = np.zeros([params.Nsnap]) # Resonant phase?
Kz = np.zeros([params.Nsnap]) # Number of buckets for spatial discretization, I think

funcs.computeUndulatorField(res_phase,Kz) # Compute undulator field

## FEL Parameters

rho1D = (1/params.gamma0)*pow( (1/8)*(params.I/params.IA)*(pow(params.K,2)/(pow(params.sigmax*params.ku,2))),(1/3))
Lgain = params.lambdau/(4*np.sqrt(3)*np.pi*rho1D)
Lsat =   params.lambdau/rho1D
Psat = 1.6*rho1D*params.Ee*params.I