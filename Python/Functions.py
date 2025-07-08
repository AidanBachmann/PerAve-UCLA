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

def primes(num): #  Sieve of Eratosthenes, adapted from: https://www.geeksforgeeks.org/python/python-program-for-sieve-of-eratosthenes/
    prime = [True for i in range(num+1)]
    pArr = []
    # Boolean array
    p = 2
    while (p * p <= num):
        # If prime[p] is not
        # changed, then it is a prime
        if (prime[p] == True):
            # Updating all multiples of p
            for i in range(p * p, num+1, p):
                prime[i] = False
        p += 1

    # Save all prime numbers
    for p in range(2, num+1):
        if prime[p]:
            pArr.append(p)
    return np.asarray(pArr)

def hammersley(D,N):
# HAMMERSLEY - Hammersley quasi-random sequence
# 
#   S = HAMMERSLEY(D,N) Generates N numbers from a D-dimensional
#   Hammersley quasirandom sequence using successive primes
#   as bases except for the last dimension which has the form
#   [1:N]'/N-1/2/N (where the last term modifies usual Hammersley
#   sequence to produce sequence in open interval (0,1)). The 
#   matrix S is of size DxN. D and N must be int.

# Copyright (c) 2008 Aki Vehtari

# This software is distributed under the GNU General Public 
# Licence (version 2 or later) please refer to the file 
# Licence.txt, included with the software, for details.
    S = np.zeros([N,D])
    S[:,-1] = (np.linspace(1,N,N,dtype='int')/N)-1/(2*N)*np.ones([N])
    pn = 2*D
    p = primes(pn)
    while(len(p) < D-1):
        pn = 2*pn
        p = primes(pn)
    P = p[0:D-1]
    for k in np.linspace(0,D-2,D-1,dtype='int'):
        pk = P[k]
        for j in np.linspace(1,N,N,dtype='int'):
            bj = j
            n = max(1,round(np.log2(bj+1)/np.log2(pk)))
            while pow(pk,n) <= bj:
                n += 1
            b = np.zeros([1,n])
            b[0,-1] = np.remainder(bj,pk)
            while (bj > 1) and (n > 1):
                n -= 1
                bj = math.floor(bj/pk)
                b[0,n-1] = np.remainder(bj,pk)
            S[j-1,k] = np.sum(np.fliplr(b)/np.power(pk,np.linspace(1,b.shape[1],b.shape[1],dtype='int')))
    return S.T

def peraveCore(oldfield,firstpass): # Push particle
    Np = params.Np # Grab number of particles
    nbins = 32 # Binning for particles?
    mpart = int(Np/nbins)
    
    radfield = np.ones([params.Nsnap,params.nslices])*params.E0 # Radiation field
    radfield[0,:] = params.profile_l*params.E0

    if firstpass == False:
        radfield[1,:] = oldfield

    thetap = np.zeros([params.Nsnap,params.nslices,Np]) # Phase space arrays
    gammap = np.zeros([params.Nsnap,params.nslices,Np])

    for islice in np.linspace(0,params.nslices-1,params.nslices,dtype='int'):
        X0 = hammersley(int(2),Np)
        gammap[0,islice,:] = params.gamma0 + params.deltagamma*X0[0,:]
        auxtheta1 = hammersley(1,mpart).T*(2*np.pi)/nbins - np.pi

        for jbin in np.linspace(0,nbins,nbins-1,dtype='int'):
            for ipart in np.linspace(0,mpart,mpart-1,dtype='int'):
                thetap[0,islice,ipart+jbin*mpart] = auxtheta1[ipart] + 2*jbin*(np.pi/nbins)
        
        if params.shotnoise > 0: # Add noise to shot
            an = 2*np.sqrt(-np.log(np.random.rand())/params.n_electron)    
            phin = np.random.rand()*2*np.pi
            for ipart in np.linspace(0,Np-1,Np,dtpye='int'):
                thetap[0,islice,ipart] -= an*np.sin(thetap[0,islice,ipart]+phin)

        '''
        if (param.prebunching ==1 )
            thetap(1,islice,:) = thetap(1,islice,:)-2.*param.bunch*sin(thetap(1,islice,:)+param.bunchphase)
        end
        if (param.prebunching < 0)
            thetab  = squeeze(thetap(1,islice,:))
                gammab = squeeze(gammap(1,islice,:))
                [thetab,gammab] = buncher(thetab,gammab,param.buncherAmp)
                for ipart = 1:Np
                thetap(1,islice,ipart) = thetab(ipart) + param.bunchphase
                gammap(1,islice,ipart) = gammab(ipart)
                end
        end

        bunching(islice) = (sum(exp(1i.*thetap(1,islice,:))/Np))
    end'''

def peravePostprocessing(): # Postprocess data from core

    return 0

def oscLoop(npasses): # Oscillator loop
    firstpass = True # Flag to indicate first pass of oscillator
    print(f'\nStarting oscillator simulation at {datetime.now().strftime("%d/%m/%Y %H:%M:%S")}.\n')

    oldfield = np.zeros([params.nslices]) # Array to store field from previous pass

    simStart = time.time() # Start time
    for i in np.linspace(0,npasses-1,npasses,dtype='int'):
        print(f'Loop {i+1} starting at {datetime.now().strftime("%d/%m/%Y %H:%M:%S")}.')
        peraveCore(oldfield,firstpass)
        firstpass = False
    simEnd = time.time() # End time