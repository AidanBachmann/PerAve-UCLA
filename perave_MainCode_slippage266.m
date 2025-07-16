%%%% PBPL FEL simulation code %%%%
%%% Input deck compatible with wafFEL, 1D period average, and GENESIS %%%
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
firstpass =1;
for jscan = 1:1
    clear power radfield thetap gammap bunch
    param.sigma_t = 40e-15;
    t1 = tic;
    tapering_strength = jscan;

Perave_User_Input_slippage266;
%% Compute the undulator field
compute_undulator_field_v5h
%% Calculate 1-D FEL parameters
rho1D = 1/param.gamma0*(1/8*param.I/IA*param.K.^2/param.sigmax^2/param.ku^2)^(1/3);
Lgain = param.lambdau/(4*sqrt(3)*pi*rho1D);
Lsat =   param.lambdau/rho1D;
Psat = 1.6*rho1D*param.Ee*param.I;
if param.tapering
    [psi1, psi2, bucket_height, capture_fraction, bucket_area, bunching_factor] = bucket_parameters(param.psir);
    a1 = 2*param.lambda0/param.lambdau*e0*param.E0/me/c^2*sin(param.psir);
    a2 = ((2*param.lambda0/param.lambdau)^1.5)*Z0*e0*param.I*sin(param.psir)^2*capture_fraction*bunching_factor/2/param.A_e/me/c^2;
    pmax_prediction=P0+param.K*(a1*lwig+a2*lwig^2/2)/(1+param.K^2)*param.Ee*param.I*capture_fraction;
    etamax(jscan) = param.K*(a1*lwig+a2*lwig^2/2)/(1+param.K^2)*capture_fraction;
    bunchlength_rms(jscan) = param.sigma_t;
    peakcurrent(jscan) = param.I;
end
%% Run the main integration routine
    
    t0 = tic;      
    perave_core_v6;
    disp(['Simulation time = ',num2str(toc(t0)./60),' min'])
    perave_postprocessor_v6   
    rad_vs_und(:,jscan) = sum(power,2)*param.lambda0*param.zsep/c;
    rad_vs_beam(:,jscan) = power(end,:);
    Eff(jscan) = Efficiency;
    if(param.itdp)
    PL(jscan) = pulselength;
    end
%
end
%% Post-process stuff
figure(100)
plot(bunchlength_rms,Eff,'b')
hold on
plot(bunchlength_rms,etamax,'r')
figure(101)
contour(rad_vs_und)
figure(102)
plot(bunchlength_rms,PL,'r')
figure(200)
plot(peakcurrent,Eff,'b')


