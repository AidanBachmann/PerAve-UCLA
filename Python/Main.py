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