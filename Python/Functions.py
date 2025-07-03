# Description: Function definitions for simulation.

# ---------- Imports ----------

import numpy as np
import math
import time
from datetime import datetime
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

def correction3D(Lgain): # Compute 3D correction to Lgain
    eta_diffr = Lgain/( (4*np.pi*pow(params.sigmax,2))/params.lambda0 )
    eta_em = (Lgain/params.betax)*((4*np.pi*params.emitx)/(params.lambda0*params.gamma0)) 
    eta_es = Lgain*(4*np.pi/params.lambdau)*params.deltagammarel

    a = [0.45, 0.57, 0.55, 1.6, 3, 2, 0.5, 2.9, 2.4, 51, 0.95, 3, 5.4, 0.7, 1.9, 1140, 2.2, 2.9, 3.2]

    eta = a[0]*pow(eta_diffr,a[1]) + a[2]*pow(eta_em,a[3]) + a[4]*pow(eta_es,a[5])
    eta += a[6]*pow(eta_em,a[7])*pow(eta_es,a[8]) + a[9]*pow(eta_diffr,a[10])*pow(eta_es,a[11]) + a[12]*pow(eta_diffr,a[13])*pow(eta_em,a[14])
    eta += a[15]*pow(eta_diffr,a[16])*pow(eta_em,a[17])*pow(eta_es,a[18])
    Lgain3D = Lgain*(1+eta)

    return Lgain3D

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

def peraveCore(): # Push particle
    Np = params.Np # Grab number of particles
    nbins = 32 # Binning for particles?
    mpart = Np/nbins 
    n_electron = (params.I*params.lambda0*params.zsep)/(params.e0*params.c) # Number of electrons?
    p1 = np.zeros([Np,1])

    tslice = np.linspace(1,params.nslices,params.nslices,dtype='int')*( (params.lambda0*params.zsep)/params.c ) # Time slices

    if (params.beamdistribution == 1): # Compute beam distribution
        profile_b = np.exp(-pow(tslice-tslice[-1]*np.ones(params.nslices)/2,2)/pow(2*params.sigma_t,2))
    else:
        profile_b = np.zeros([params.nslices])
        profile_b[abs(tslice-tslice[-1]*np.ones([params.nslices])/2)<params.sigma_t] = 1

    if (params.laserdistribution == 1): # Compute laser distribution
        profile_l = np.exp(-pow(tslice-params.slippage*np.ones([params.nslices]),2)/(2*pow(params.sigma_l,2)))
    else:
        profile_l = np.zeros([params.nslices])
        profile_l[abs(tslice-params.slippage*np.ones([params.nslices]))<params.sigma_l] = 1
    

def peravePostprocessing(): # Postprocess data from core

    return 0

def oscLoop(npasses): # Oscillator loop
    print(f'\nStarting oscillator simulation at {datetime.now().strftime("%d/%m/%Y %H:%M:%S")}.\n')
    simStart = time.time() # Start time
    for i in np.linspace(0,npasses-1,npasses,dtype='int'):
        print(f'Loop {i+1} starting at {datetime.now().strftime("%d/%m/%Y %H:%M:%S")}.')
        peraveCore()
    simEnd = time.time() # End time