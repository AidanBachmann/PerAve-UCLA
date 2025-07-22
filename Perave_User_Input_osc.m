% Perave_code_user input
%%%%User entered parameters%%%%%
%% from Duris et al. TESSO paper 

% Set random number generator seed
rng(3141592653)

%% Undulator parameters
param.lambdau = 3.14e-2;                                 % undulator period (m)
param.K = 1.934; %e0*Bfield*/me/c/ku;               % RMS undulator parameter
param.ku = 2.*pi./param.lambdau;                   % undulator wavenumber
lwig = 9.42;                                                             % Undulator length (m), paper defines number of undulator periods and lambdau, so lwig is computed from these vals 
param.undulator_type = 0; % Type of undulator (0 for planar, 1 for helical)
% Tapering options
param.tapering = 0;                                         % tapering (-1 acceleration ; 0 no tapering ; 1 decelation)    
param.z0 = 0;
param.psir = pi/6;

%% Simulation control options
param.suppress_plots = 1; % Suppress plotting in perave_postprocessor_v6 (1 for True, 0 for False)
param.phasespacemovie = 0;
param.itdp = 0; % 1 for time-dependent simulation, 0 for time independent
param.prebunching = 1;                                                                  % set to 1 to start from a pre-bunched beam. 
param.changeresphase = 0;
saveoutput=0;
% Set simulation length and # of snapshots
param.delz=1;
param.und_periods = round(lwig/param.lambdau);                         % number of undulator periods to simulate
param.Nsnap = round(lwig/param.lambdau/param.delz);                % number of snapshots to take over the length of the undulator
param.zsep = 5;                                                              
Nslip=round(param.Nsnap/param.zsep);
param.shotnoise = 1;
param.lambda0 = 13.5e-9;                                    % seed wavelength (m)
param.k = 2*pi/param.lambda0;                                     % wavenumber in free space
param.nslices = 4*Nslip+4*round(param.sigma_t/param.zsep/param.lambda0*c);

if(~param.itdp)
    param.nslices = 1;
end
if(param.itdp)
    param.Nsnap = floor(param.und_periods/param.delz);        % Note that the result must be an integer
end
param.stepsize = param.lambdau*param.delz;

%% Electron beam parameters
gamma0 = sqrt(param.k/2/param.ku*(1+param.K^2));param.gamma0=gamma0;          % relativistic gamma factor
param.Np = 512;                                          % # of macroparticles (500-1000 well)
param.Ee = gamma0*me*c^2/e0;                  % Total e-beam energy (eV)
energyspread = 1*20e-15/param.sigma_t;                                       % Absolute energy spread MeV 
param.deltagammarel = energyspread/gamma0/0.511;          % Relative energy spread dgamma/gamma should be 10^-4
param.deltagamma = gamma0*param.deltagammarel;
param.bunch = 0.000001;                                  % Initial bunching factor
param.bunchphase = -param.psir-pi/2;                     % Initial bunching phase
param.buncherAmp = 5;

betax = 5.1; % Average beta function, 4.7544975 (something weird going on here)
emitx = 1e-6; % normalized emittance in m*rad
charge = param.charge; % Computed from peak current and sigma_t
%param.sigma_t = 40e-15;
if (param.beamdistribution == 1)
    param.I = charge/sqrt(2*pi)/param.sigma_t;     % beam current 
else
    param.I = charge/2/param.sigma_t;              % beam current 
end

param.sigmax = sqrt(betax*emitx/gamma0);            % beam radius
param.A_e = 2*pi*param.sigmax^2;                          % beam cross section 
Simulation_temporal_window=param.nslices*param.zsep*param.lambda0/c;

%% radiation parameters
P0 = 1e3; param.P0=P0;                                               % Peak input power (W) 
A_mode = param.A_e;                                                     % 1D code. area is same for e_beam and radiation
param.waist = sqrt(A_mode*2/pi);
zr = pi*param.waist^2/param.lambda0;                          % Rayleigh length of seed (m)
param.E0 = sqrt(2*P0/c/eps0/A_mode/2);                        % Assume circular polarization
param.slippage = param.nslices/2*param.lambda0*param.zsep/c;
param.sigma_l = 2400e-15;

%% Simplifying constants
param.chi2 = e0/me/c^2;
param.chi1 = mu0*c/2*param.I/param.A_e;
param.chi = (param.K.^2)/(2*(1 + param.K.^2)); % Argument of the JJ factort that appears in the coupling factor (this changes with tapering)
param.fc = besselj(0,param.chi) - besselj(1,param.chi); % Coupling constant

%% 1D FEL parameters
param.omega0 = param.ku*c;
if param.undulator_type == 0
    param.rho1D = (1/param.gamma0*(1/8*param.I/IA*param.K.^2/param.sigmax^2/param.ku^2)^(1/3))*param.fc.^(2/3); % 1D FEL parameter
else
    param.rho1D = 1/param.gamma0*(1/8*param.I/IA*param.K.^2/param.sigmax^2/param.ku^2)^(1/3); % 1D FEL parameter
end
param.Lgain = param.lambdau/(4*sqrt(3)*pi*param.rho1D); % Gain length
param.Lsat =   param.lambdau/param.rho1D; % Saturation length
param.Psat = 1.6*param.rho1D*param.Ee*param.I; % Saturation power
param.gammar = sqrt( param.k*(1 + param.K.^2)/(2*param.ku) ); % Resonant energy
param.sigma = (4*param.rho1D.^2)*(1+param.K.^2)/param.K.^2; % Sigma parameter
param.delta = (param.gamma0.^2 - param.gammar.^2)/(2*(param.gammar.^2)*param.rho1D); % Delta parameter
param.lambda = findRoots(param.rho1D,param.sigma,param.delta); % Exponential gain factor, used to compute gain length
param.gain_lambda = param.lambdau/(param.lambda*4*pi*param.rho1D);