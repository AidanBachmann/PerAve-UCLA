# Description: Function definitions for simulation.

# ---------- Imports ----------

import numpy as np
import math
import matplotlib.pyplot as plt
import Params as params

# ---------- Function Definition ----------

def computeUndulatorField(res_phase,Kz): # Compute undulator field, I guess
    if params.tapering != 0: # Check if using tapering
        startindex = max(0,math.floor(params.z0/params.stepsize)) # Discretized initial position 
        res_phase[startindex:params.Nsnap] = params.psir # Set phase
    Kz[0] = params.K # Store initial RMS K value in Kz array

def bucket_parameters(): # Compute bucket parameters from params.psir. Definition coming later.
    bucket_parameters = np.zeros([7]) # Initialize array to store bucket parameters
    
    return bucket_parameters

def filter(sigma_omega,plotFilter): # Compute filters
    jfreq = np.linspace(0,params.nslices,params.nslices,dtype='int')
    filter = np.exp(-pow(jfreq-params.nslices*np.ones(params.nslices)/2,2)/(2*pow(sigma_omega,2)))
    filter2 = np.zeros([len(jfreq)],dtype='complex') # Complex transfer function
    filter3 = np.zeros([len(jfreq)],dtype='complex')

    for jf in np.linspace(0,params.nslices-1,params.nslices,dtype='int'):
        y = ((jf+1) - params.nslices/2)/sigma_omega
        if(y >= 1):
            filter2[jf] = y - np.sqrt(pow(y,2) - 1)
        elif(y<=-1):
            filter2[jf] = y + np.sqrt(pow(y,2) - 1)
        else:
            filter2[jf] = y + 1j*np.sqrt(1 - pow(y,2))
        omega_m = params.nslices/2
        Q = 1/sigma_omega
        filter3[jf] = ( 1j*((jf+1)/Q) )/( pow(omega_m,2) - pow((jf+1),2) + 1j*((jf+1)/Q) )
    filterdelay = round(params.nslices/(2*np.pi*sigma_omega))
    #print(f'Filter = {filter}\nFilter2 = {filter2}\nFilter3 = {filter3}\nFilter Delay = {filterdelay}\n')
    if plotFilter:
        plt.scatter(jfreq,np.real(filter3),label='Re[filter3]')
        plt.scatter(jfreq,np.imag(filter3),label='Im[filter3]')
        plt.xlabel('Frequencies')
        plt.ylabel('Filter')
        plt.legend()
        plt.show()
    return filter3