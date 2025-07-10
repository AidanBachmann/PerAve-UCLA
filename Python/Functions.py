# Description: Function definitions for simulation.

# ---------- Imports ----------

import numpy as np
import scipy
import math
import time
from datetime import datetime
from numba import jit
import matplotlib.pyplot as plt
import Params as params

# ---------- Function Definition ----------

def printParams(Lgain3D): # Print parameters to terminal (unfortunately, hardcoded)
    print('Simulation Parameters\n')
    print(f'sigma_t: {params.sigma_t}\nuse3DCorrection: {params.use3Dcorrection}\nbeamdistribution: {params.beamdistribution}\nlaserdistribution: {params.laserdistribution}')
    print(f'lambdau: {params.lambdau}\nK: {params.K}\nku: {params.ku}\ntapering: {params.tapering}\nz0: {params.z0}\npsir: {params.psir}\nphasespacemovie: {params.phasespacemovie}')
    print(f'itdp: {params.itdp}\nprebunching: {params.prebunching}\nchangeresphase: {params.changeresphase}\ndelz: {params.delz}\nund_periods: {params.und_periods}')
    print(f'Nsnap: {params.Nsnap}\nzsep: {params.zsep}\nshotnoise: {params.shotnoise}\nlambda0: {params.lambda0}\nk: {params.k}\nnslices: {params.nslices}')
    print(f'stepsize: {params.stepsize}\ngamma0: {params.gamma0}\nNp: {params.Np}\nEe: {params.Ee}\ndeltagammarel: {params.deltagammarel}\ndeltagamma: {params.deltagamma}')
    print(f'bunch: {params.bunch}\nbunchphase: {params.bunchphase}\nbuncherAmp: {params.buncherAmp}\nI: {params.I}\nLgain3D: {Lgain3D}')

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

@jit(nopython = True)
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

@jit(nopython = True)
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
    S = np.zeros((N,D))
    S[:,-1] = (np.linspace(1,N,N).astype('int')/N)-1/(2*N)*np.ones((N))
    pn = 2*D
    p = primes(pn)
    while(len(p) < D-1):
        pn = 2*pn
        p = primes(pn)
    P = p[0:D-1]
    for k in np.linspace(0,D-2,D-1).astype('int'):
        pk = P[k]
        for j in np.linspace(1,N,N).astype('int'):
            bj = j
            n = max(1,round(np.log2(bj+1)/np.log2(pk)))
            while pow(pk,n) <= bj:
                n += 1
            b = np.zeros((1,n))
            b[0,-1] = np.remainder(bj,pk)
            while (bj > 1) and (n > 1):
                n -= 1
                bj = math.floor(bj/pk)
                b[0,n-1] = np.remainder(bj,pk)
            S[j-1,k] = np.sum(np.fliplr(b)/np.power(pk,np.linspace(1,b.shape[1],b.shape[1]).astype('int')))
    return S.T

@jit(nopython = True)
def buncher(thetab,gammab,amp): # Buncher
    amp1 = amp/4.7
    BR56 = np.pi/(2*(amp+3))
    BR561 = np.pi/(amp1+3)

    #plt.figure(3)
    gamma0 = np.mean(gammab)
    gammaspread = np.std(gammab)
    #plt.plot(thetab,gammab)
    gammarel = (gammab - gamma0)/gammaspread
    gammarel = gammarel - amp1*np.sin(thetab)
    #plt.figure(6)
    #plt.plot(thetab,gammarel)
    phaseb = thetab + gammarel*BR561
    #plt.figure(4)
    #plt.plot(phaseb,gammarel)
    gammarel = gammarel - amp*np.sin(phaseb)
    phaseb = phaseb + gammarel*BR56
    energyb = gammarel*gammaspread + gamma0
    #plt.figure(5)
    #plt.plot(phaseb,energyb)
    return phaseb,energyb

@jit(nopython = True)
def push_FEL_particles_RK4(phasespace,evalue,kvalue,chi1):
    gammar_sq = params.lambdau/(2*params.lambda0)*(1+pow(kvalue,2))
    sc = 1
    # Euler method for field (for some reason the most accurate...)
    newevalue = evalue - params.stepsize*(chi1*kvalue*sc*np.mean(np.exp(-1j*phasespace[:,0])/phasespace[:,1]))

    # Leapfrog method for the field

    # RK-2 for the field
    #k1e=-1*chi1*kvalue*mean(exp(-1j*phasespace[:,0])./phasespace[:,1])
    #y1e=evalue+k1e*param.stepsize/2
    #k2e=-1*chi1*kvalue*mean(exp(-1j*(phasespace[:,0]+param.stepsize/2))./(phasespace[:,1]+param.stepsize/2))
    #newevalue=evalue+k2e*param.stepsize
    
    # RK-4 for the particles

    k1theta = params.stepsize*(params.ku*(1-(gammar_sq/np.square(phasespace[:,1]))))
    k1gamma = params.stepsize*(params.chi2*(kvalue/phasespace[:,1])*np.real(evalue*sc*np.exp(1j*phasespace[:,0])))
    
    k2theta = params.stepsize*(params.ku*(1-(gammar_sq/np.square(phasespace[:,1]+0.5*k1gamma))))
    k2gamma = params.stepsize*(params.chi2*(kvalue/(phasespace[:,1]+0.5*k1gamma))*np.real(evalue*sc*np.exp(1j*(phasespace[:,0]+0.5*k1theta))))
    
    k3theta = params.stepsize*(params.ku*(1-(gammar_sq/np.square(phasespace[:,1]+0.5*k2gamma))))
    k3gamma = params.stepsize*(params.chi2*(kvalue/(phasespace[:,1]+0.5*k2gamma))*np.real(evalue*sc*np.exp(1j*(phasespace[:,0]+0.5*k2theta))))
    
    k4theta = params.stepsize*(params.ku*(1-(gammar_sq/np.square(phasespace[:,1]+k3gamma))))
    k4gamma = params.stepsize*(params.chi2*(kvalue/(phasespace[:,1]+k3gamma))*np.real(evalue*sc*np.exp(1j*(phasespace[:,0]+k3theta))))
    
    newphasespace = np.vstack((phasespace[:,0] + 1/6*(k1theta+2*k2theta+2*k3theta+k4theta),phasespace[:,1] + 1/6*(k1gamma+2*k2gamma+2*k3gamma+k4gamma))).T

    return newphasespace,newevalue

@jit(nopython = True)
def bukh(phi):
    return np.sqrt(np.cos(phi)-(np.pi/2-phi)*np.sin(phi))

@jit(nopython = True)
def peraveCore(oldfield,firstpass,Kz,res_phase): # Push particle
    Np = params.Np # Grab number of particles
    nbins = 32 # Binning for particles?
    mpart = int(Np/nbins)
    
    radfield = np.ones((params.Nsnap,params.nslices)).astype('complex')*params.E0 # Radiation field
    radfield[0,:] = params.profile_l*params.E0

    if firstpass == False:
        radfield[0,:] = oldfield

    thetap = np.zeros((params.Nsnap,params.nslices,Np)) # Phase space arrays
    gammap = np.zeros((params.Nsnap,params.nslices,Np))
    bunching = np.zeros((params.nslices)).astype('complex')

    for islice in np.linspace(0,params.nslices-1,params.nslices).astype('int'):
        X0 = hammersley(int(2),Np)
        gammap[0,islice,:] = params.gamma0 + params.deltagamma*X0[0,:]
        auxtheta1 = hammersley(1,mpart).T*(2*np.pi)/nbins - np.pi

        for jbin in np.linspace(0,nbins-1,nbins).astype('int'):
            for ipart in np.linspace(0,mpart-1,mpart).astype('int'):
                thetap[0,islice,ipart+jbin*mpart] = auxtheta1[ipart,0] + 2*jbin*(np.pi/nbins)
        
        if params.shotnoise > 0: # Add noise to shot
            an = 2*np.sqrt(-np.log(np.random.rand())/params.n_electron)    
            phin = np.random.rand()*2*np.pi
            for ipart in np.linspace(0,Np-1,Np).astype('int'):
                thetap[0,islice,ipart] -= an*np.sin(thetap[0,islice,ipart]+phin)

        if params.prebunching == 1:
            thetap[0,islice,:] -= 2*params.bunch*np.sin(thetap[0,islice,:] + params.bunchphase)
        
        if params.prebunching < 0:
            thetab = thetap[0,islice,:]
            gammab = gammap[0,islice,:]
            thetab,gammab = buncher(thetab,gammab,params.buncherAmp)
            for ipart in np.linspace(0,Np-1,Np).astype('int'):
                thetap[0,islice,ipart] = thetab[ipart] + params.bunchphase
                gammap[0,islice,ipart] = gammab[ipart]

        bunching[islice] = np.sum(np.exp(1j*thetap[0,islice,:])/Np)
    
    print(f'Finished init.')

    ## Solve system of equations
    res_step = params.und_periods*params.lambdau/params.Nsnap   
    total_simtime = 0
    hl = 0
    z = np.zeros((1))
    gammares = np.zeros((params.Nsnap))
    gammares[0] = np.sqrt(params.lambdau*(1+pow(Kz[0],2))/(2*params.lambda0)) 
    # Constant for the resonant phase based tapering   
    const_resp = (1/params.chi2)*(params.lambdau/(2*params.lambda0))
    slip = 0
    bunch = np.zeros((params.Nsnap,params.nslices)).astype('complex') # Bunching array

    if params.itdp == int(1): # Time dependent simulation
        for ij in np.linspace(0,params.Nsnap-2,params.Nsnap-1).astype('int'):  # Takes Nsnap snapshots along length of undulator
            for islice in np.linspace(0,params.nslices-1,params.nslices).astype('int'):
                gammaf = gammap[ij,islice,:]
                thetaf = thetap[ij,islice,:]
                E_q0 = radfield[ij,islice]
                chi1 = (params.mu0*params.c/2)*(params.I*params.profile_b[islice]/params.A_e)
                # RK4th order integration     
                phasespaceold = np.vstack((thetaf,gammaf)).T
                evaluesold = E_q0
                
                phasespacenew,evaluesnew = push_FEL_particles_RK4(phasespaceold,evaluesold,Kz[ij],chi1)
                thetap[ij+1,islice,:] = phasespacenew[:,0]
                gammap[ij+1,islice,:] = phasespacenew[:,1]
                radfield[ij+1,islice] = evaluesnew
            # Slippage of the radiation field 
            if ( np.mod(ij,params.zsep) == (params.zsep-1) ): # ***** Check on indexing
                B = radfield[ij+1,:]
                B = np.roll(B,1)
                radfield[ij+1,:] = B
                if (firstpass == False):
                    radfield[ij+1,0] = 0
                else:
                    radfield[ij+1,0] = params.E0*params.profile_l[0]
            
            for k in np.linspace(0,params.nslices-1,params.nslices).astype('int'): # Loop to compute mean since jit doesn't support axis kwarg in mean function
                bunch[ij,k] = np.mean(np.exp(1j*thetap[ij,k,:])) # shape(thetap) = (Nsnap,nslices,Np)
            
            if firstpass == True:
                if params.tapering_strength == 0:
                    Klz = max(np.abs(radfield[0,:]))
                elif params.tapering_strength == 1:
                    Klz = max(np.abs(radfield[ij,:]))
                elif params.tapering_strength == 2:
                    Klz = np.mean(np.abs(radfield[ij,:]))
                Kz[ij+1] = Kz[ij] - (params.stepsize/const_resp)*Klz*np.sin(res_phase[ij]) # ***** I think this should be outside of the firstpass if statement. Pietro code has it inside if statement.
            
            gammares[ij+1] = np.sqrt(params.lambdau*(1+pow(Kz[ij],2))/(2*params.lambda0))
            print(f'Finished snapshot {ij+1} out of {params.Nsnap-1}')
        # Remove slices within one total slippage length
        '''radfield[:,0:params.Nslip-1] = [] # Fix later
        gammap[:,0:params.Nslip-1,:] = []
        thetap[:,0:params.Nslip-1,:] = []
        params.profile_l[0:params.Nslip-1] = []
        params.profile_b[0:params.Nslip-1] = []
        bunch[:,0:params.Nslip-1] = []'''
    else: # Time independent simulation
        deltagammamax = 1
        for ij in np.linspace(0,params.Nsnap-2,params.Nsnap-1).astype('int'):  # Takes Nsnap snapshots along length of undulator
            gammaf = gammap[ij,0,:]
            thetaf = thetap[ij,0,:]
            E_q0 = radfield[ij,0]
            
            # RK4th order integration     
            phasespaceold = np.vstack((thetaf,gammaf)).T
            evaluesold = E_q0
            phasespacenew,evaluesnew = push_FEL_particles_RK4(phasespaceold,evaluesold,Kz[ij],chi1)     
            thetap[ij+1,0,:] = phasespacenew[:,0]
            gammap[ij+1,0,:] = phasespacenew[:,1]
            radfield[ij+1,0] = evaluesnew
            
            # Compute undulator field at next step (constant res phase)
            Kz[ij+1] = Kz[ij] - (params.stepsize/const_resp)*np.abs(radfield[ij,0])*np.sin(res_phase[ij])
            gammares[ij+1] = np.sqrt(params.lambdau*(1 + pow(Kz[ij],2))/(2*params.lambda0))
            bunch[ij] = np.mean(np.exp(1j*thetap[ij,0,:]))

            '''if (ij > 40000) and (params.changeresphase):
                deltagamma = np.sqrt(Kz[ij+1]*np.mean(np.abs(radfield[ij,:])))*bukh(res_phase[ij])
                deltagammamax = max(deltagammamax,deltagamma)
                if ( deltagamma < deltagammamax):
                     #newphase = fsolve(@(phi) sqrt(Kz(ij+1).*mean(abs(radfield(ij,:)),2))*bukh(phi)-deltagammamax, res_phase(ij)); # Not sure I can make this compatible with jit, I'll leave it commented out for now
                     pass'''
            print(f'Finished snapshot {ij+1} out of {params.Nsnap-1}')
    print('Finished solving eqs.')

    # Calculate radiation power 
    power = np.power(np.abs(radfield),2)/(params.Z0*params.A_e)

    return power,radfield,gammap,thetap

def spectrum_calc(field,xlamds,zsep): # Compute spectrum
    omegas = (2*np.pi*params.c)/xlamds
    df = (2*np.pi*params.c)/(params.nslices*zsep*xlamds)
    ft = scipy.fft.fftshift(scipy.fft.fft(field))
    omega = df*np.linspace(1,len(ft),len(ft),dtype='int')
    omega -= np.median(omega) + omegas
    omega /= omegas
    power_spectrum = np.power(abs(ft),2)
    omega -= 1
    return power_spectrum,omega

def peravePostprocessing(radfield,power,gammap,thetap,rho1D): # Postprocess data from core
    kw = 2*np.pi/params.lambdau
    zpos = np.linspace(0,params.Nsnap-1,params.Nsnap,dtype='int')*params.stepsize # ***** Matlab 1 indexing means that Matlab version runs over np.linspace(1,params.Nsnap,params.Nsnap,dtype='int')

    ## Spectrum as a function of z
    fundpower = np.zeros([params.Nsnap]) # Power
    sidebandpower = np.zeros([params.Nsnap])
    if params.itdp == 1:
        omegamin = -10e-4
        omegamax = 10e-4
        fig,ax = plt.subplots(nrows=1,ncols=2,figsize=(16,9))
        ax = ax.flatten()
        fig.suptitle('Power Spectrum')
        xax = (params.zsep*params.lambda0*1e15/params.c)*np.linspace(0,power.shape[1]-1,power.shape[1],dtype='int')
        for n in np.linspace(0,params.Nsnap-1,params.Nsnap,dtype='int'):
            powerspec,omega = spectrum_calc(radfield[n,:],params.lambda0,params.zsep)
            idxMin = np.where(omega > omegamin) # Sideband index
            idxMax = np.where(omega < omegamax)
            try:
                fundspectrum = powerspec[idxMin[0][0]:idxMax[0][-1]]
                fundpower[n] = np.trapz(fundspectrum)/np.trapz(powerspec)
            except:
                pass
        
        ax[0].plot(omega/rho1D,abs(powerspec))
        ax[0].set_xlabel(r'$\frac{\delta\omega}{\rho\omega}$')
        ax[0].set_ylabel(r'P($\omega$) [arb. units]')    
        ax[0].set_yscale('log')
        ax[0].set_xlim([-20,20])
    
        ax[1].plot(xax,power[-1,:])
        ax[1].set_xlim([0,max(xax)*1.025])
        ax[1].set_xlabel('t [fs]')
        ax[1].set_ylabel('Output Radiation Power [W]')

        plt.show()
        
    
    ## Radiation Power and spectrum at exit
    fig,ax = plt.subplots(nrows=2,ncols=3,figsize=(16,9))
    ax = ax.flatten()
    fig.suptitle('Simulation Output')

    ax[0].plot(zpos,np.mean(power,axis=1),color='b',label='Avg')
    ax[0].plot(zpos,np.max(power,axis=1),color='r',label='Max')
    ax[0].set_xlim([0,zpos[-1]*1.025])
    ax[0].set_yscale('log')
    ax[0].set_xlabel('z Position')
    ax[0].set_ylabel('Power')
    ax[0].set_title('Radiation Power along the beam')
    ax[0].legend()

    if params.itdp == 1:
        xax = (params.zsep*params.lambda0*1e15/params.c)*np.linspace(0,power.shape[1]-1,power.shape[1],dtype='int')
        ax[1].plot(xax,power[-1,:],label='Final Power')
        ax[1].plot(xax,(np.max(power[-1,:]))*params.profile_l,linestyle='-',label=f'Norm Initial {params.P0/1e9} GW')
        ax[1].plot(xax,(np.max(power[-1,:])/2)*params.profile_b,marker='.',label=f'Current Profile {params.I/1e3} kA')
        ax[1].set_xlim([0,np.max(xax)*1.025])
        ax[1].set_xlabel('t [fs]')
        ax[1].set_ylabel('Power [W]')
        ax[1].set_title('Power as a function of time')
        ax[1].legend()

        powerspec,omega = spectrum_calc(radfield[-1,:],params.lambda0,params.zsep)
        xax = (omega+1)*params.hbar*2*np.pi*params.c/params.lambda0
        ax[2].plot(xax,powerspec,color='b')
        ax[2].set_xlim(min(xax),max(xax))
        ax[2].set_yscale('log')
        ax[2].set_xlabel('Photon Energy [eV]')
        ax[2].set_ylabel(r'P($\omega$) [arb. units]')
        ax[2].set_title('Output Spectrum')

    meanenergy = np.zeros([params.Nsnap])
    for ij in np.linspace(0,params.Nsnap-1,params.Nsnap,dtype='int'):
        meanenergy[ij] = sum(np.mean(gammap[ij,:,:],axis=2)*params.profile_b)/sum(params.profile_b)
    
    xax = params.stepsize*np.linspace(0,params.Nsnap-1,params.Nsnap,dtype='int')
    ax[3].plot(xax,meanenergy)
    ax[3].set_xlabel('z')
    ax[3].set_ylabel(r'$\gamma$')
    ax[3].set_xlim([0,params.Nsnap*1.025])

    plt.show()
    input('WAIT')
        


def oscLoop(npasses,Kz,res_phase,rho1D): # Oscillator loop
    firstpass = True # Flag to indicate first pass of oscillator
    print(f'\nStarting oscillator simulation at {datetime.now().strftime("%d/%m/%Y %H:%M:%S")}.\n')

    oldfield = np.zeros([params.nslices]) # Array to store field from previous pass

    simStart = time.time() # Start time
    for i in np.linspace(0,npasses-1,npasses,dtype='int'):
        print(f'Loop {i+1} starting at {datetime.now().strftime("%d/%m/%Y %H:%M:%S")}.')
        start = time.time()
        power,radfield,gammap,thetap = peraveCore(oldfield,firstpass,Kz,res_phase)
        end = time.time()
        print(f'Finished loop {i+1} in {end-start} seconds.\n')
        peravePostprocessing(radfield,power,gammap,thetap,rho1D)
        firstpass = False
    simEnd = time.time() # End time
    print(f'Finished simulation in {simEnd-simStart} seconds.\n')