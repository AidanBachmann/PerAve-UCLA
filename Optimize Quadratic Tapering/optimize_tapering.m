%% Optimize coefficients a, b, and z0 for quadratic tapering to maximize power output.
clear all
close all

%% physical constants
global c mu0 e0 me eps0 IA Z0

c = 2.99792458.*10^8;                                            % speed of light
e0 = 1.60217657e-19;                                             % electron charge
me = 9.10938291e-31;                                             % electron mass
eps0 = 8.85418782e-12;                                           % eps_0
mu0 = 1.256637e-6;                                               % mu_0
IA = 17045;                                                      % Alfven current
Z0 = 376.73;                                                     % Impedance of free space         

%% Load the User Determined initial conditions
clear power radfield thetap gammap bunch
param.sigma_t = 0.9e-12;
param.charge = 2255.965e-12; % Charge
param.use3Dcorrection  = 1;
param.beamdistribution = 1;       % Using GENESIS flag: 2-uniform 1-gaussian
param.laserdistribution = 1;         % Using GENESIS flag: 2-uniform 1-gaussian
recirculate = 0;
t1 = tic;
Perave_User_Input_osc;

x0 = [param.a0,param.b0,param.zs]; % Initial guessses for optimal tapering parameters
res = fminsearch(@sim,x0); % Optimize parameters
disp('Optimal parameters for quadratic tapering:');
disp(res);