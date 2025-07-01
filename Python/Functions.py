# Description: Function definitions for simulation.

# ---------- Imports ----------

import numpy as np
import Params as params
import math

# ---------- Function Definition ----------

def computeUndulatorField(res_phase,Kz): # Compute undulator field, I guess
    if params.tapering != 0: # Check if using tapering
        startindex = max(0,math.floor(params.z0/params.stepsize)) # Discretized initial position 
        res_phase[startindex:params.Nsnap] = params.psir # Set phase
    Kz[0] = params.K # Store initial RMS K value in Kz array

def bucket_parameters(): # Compute bucket parameters from params.psir. Definition coming later.
    bucket_parameters = np.zeros([7]) # Initialize array to store bucket parameters
    
    return bucket_parameters

def filter(sigma_omega): # Compute filters
    jfreq = np.linspace(1,params.nslices,params.nslices,dtype='int')
    filter = np.exp(-pow(jfreq-params.nslices/2,2)/(2*pow(sigma_omega,2)))
    filter2 = np.zeros([len(jfreq)],dtype='complex') # Complex transfer function

    for jf in np.linspace(0,params.nslices-1,params.nslices,dtype='int'):
        y = (jfreq-params.nslices/2)/sigma_omega

    return 0