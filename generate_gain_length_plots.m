%% Create plots for gain length as a function of current in exponential gain regime.
%% Loop over current values, compute numerical and theoretical gain lengths, plot L(I) and error.
%% Current is computed by setting charge.
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
charge_i = 0.10*2255.965e-12; % Initial charge to test
charge_f = 10.0*2255.965e-12; % Final charge to test
Nc_steps = 100; % Number of steps to take between charge_i and charge_f
chargeArr = linspace(charge_i,charge_f,Nc_steps); % Array of charge values
currArr = (1:Nc_steps)*0; % Array to store current values
LgainTarr = (1:Nc_steps)*0; % Array to store theoretical gain lengths
LgainNarr = (1:Nc_steps)*0; % Array to store numerical gain lengths
t1 = tic;
for jdx = 1:Nc_steps
    param.charge = chargeArr(jdx); % Set charge for loop
    Perave_User_Input_osc;
    currArr(jdx) = param.I; % Store current
    %% Compute the undulator field
    compute_undulator_field_v5h;
    
    %% Calculate 1-D FEL parameters
    if param.tapering
        [psi1, psi2, bucket_height, capture_fraction, bucket_area, bunching_factor] = bucket_parameters(param.psir);
        a1 = 2*param.lambda0/param.lambdau*e0*param.E0/me/c^2*sin(param.psir);
        a2 = ((2*param.lambda0/param.lambdau)^1.5)*Z0*e0*param.I*sin(param.psir)^2*capture_fraction*bunching_factor/2/param.A_e/me/c^2;
        pmax_prediction=P0+param.K*(a1*lwig+a2*lwig^2/2)/(1+param.K^2)*param.Ee*param.I*capture_fraction;
        etamax = param.K*(a1*lwig+a2*lwig^2/2)/(1+param.K^2)*capture_fraction;
        bunchlength_rms = param.sigma_t;
        peakcurrent = param.I;
    end
    calculate_3Dcorrection; 
    
    %% Run the main integration routine
    cavitydetuning = -16;    % In units of zsep
    transmission = 0.66;      % Power transmission through one cavity pass 
                                          % losses = 1 - transmission                                      
    sigma_omega = 0.003*param.nslices*param.zsep;     % Filter fractional bandwidth. 
    firstpass =1;
    tapering_strength = 2;   % 0 max of slices at time 0 
                                          % 1 max of slices
                                          % 2 avg of slices
    %% Filter definition (Filter2 is a complex transfer function. Cavity detuning needs to be adjusted to 12)
        jfreq = 1:param.nslices;
        filter = exp(-(jfreq-param.nslices/2).^2/2/sigma_omega^2);
        for jfreq = 1:param.nslices;
        y = (jfreq-param.nslices/2)/sigma_omega;
        if(y>=1)
            filter2(jfreq) = y-sqrt(y.^2-1);
        elseif(y<=-1)
            filter2(jfreq) = (y+sqrt(y.^2-1));
        else
            filter2(jfreq) = y+i*sqrt(1-y.^2);
        end
            omega_m=param.nslices/2;
            Q = 1/sigma_omega;
            filter3(jfreq) = 1i*jfreq/Q / (omega_m^2-jfreq^2+1i*jfreq/Q);   
        end    
        filterdelay = round(param.nslices/2/pi/sigma_omega);
        figure(200)
        plot(filter)
        hold on
        plot(abs(filter2))
        plot(abs(filter3),'k')
        hold off
        figure(201)
        plot(angle(filter2))
        hold on
        plot(angle(filter3),'k')
        hold off
    
     %% Oscillator loop
    for npasses = 1:1
        clear power radfield thetap gammap bunch
        t0 = tic;
        perave_core_v6;
        disp(['Simulation time = ',num2str(toc(t0)./60),' min'])
        perave_postprocessor_v6   
        rad_vs_und(:,npasses) = sum(power,2)*param.lambda0*param.zsep/c;
        rad_vs_beam(:,npasses) = power(end,:);
        Eff(npasses) = Efficiency;
        PL(npasses) = pulselength;
        oldfield(1:param.nslices) =0;
        
        if cavitydetuning>0
        oldfield(1,cavitydetuning+1:cavitydetuning+size(radfield,2)) = radfield(end,:)*sqrt(transmission);
        else
        oldfield(1,1:1+cavitydetuning+size(radfield,2)) = radfield(end,-cavitydetuning:end)*sqrt(transmission);    
        end
        pause(0.5)
    
        %%
        pltspec = abs(fftshift(fft(oldfield)).*filter3);
        filterfield = ifft(ifftshift(fftshift(fft(oldfield) ).*filter3));
        if param.suppress_plots == 0
            figure(8)
            subplot(1,2,1)
            plot(pltspec);
            subplot(1,2,2)
            plot(power(end,:),'k')
            hold on
            plot(abs(filterfield).^2/377*param.A_e,'g')
            plot(abs(oldfield).^2/377*param.A_e,'r')
            plot(profile_b*max(power(end,:))*0.5,'b')
            hold off
            pause(0.5);
        end
        oldfield = filterfield;
        firstpass = 0;                                  % Start recirculation
    end
    %% Post-process stuff
    figure(100)
    plot(max(rad_vs_und),'b')
    figure(101)
    plot([1:1:param.Nsnap]*param.stepsize,rad_vs_und(:,end),'r')
    hold on
    plot([1:1:param.Nsnap]*param.stepsize, meanenergy*charge*511000)
    xlim([0,param.Nsnap*param.stepsize])
    title('Radiation energy along undulator')
    figure(102)
    plot(PL)
    title 'pulselength'
    figure(103)
    plot(Eff)
    
    figure(300)
    try
        contourf([1:size(rad_vs_beam,1)]*param.zsep*param.lambda0/c,[1:1],rad_vs_beam');
    catch
        disp('Error using contourf, argument must be at least 2x2.')
    end
    LgainTarr(jdx) = param.Lgain; % Store theoretical gain length
    LgainNarr(jdx) = param.gainlen; % Store numerical gain length
end
normErr = abs(LgainTarr(1,:) - LgainNarr(1,:))./LgainTarr(1,:); % Compute normalized error
%% Generate plots
figure(27)
subplot(1,2,1);
plot(currArr./1000,LgainTarr,color='g'); % Plot theoretical gain length as a function of current
hold on
scatter(currArr./1000,LgainNarr,color='r'); % Plot numerical gain length as a function of current
title('Gain Length as a Function of Current');
xlim([min(currArr./1000)*0.95,max(currArr./1000)*1.05]);
ylim([min(LgainTarr)*0.975,max(LgainTarr)*1.025]);
xlabel('Current (kA)')
ylabel('Gain Length (m)')
legend('Theoretical Gain Length','Numerical Gain Length');
hold off
subplot(1,2,2);
scatter(currArr./1000,normErr);
set(gca,'yscale','log')
hold on
title('Normalized Error in Gain Length as a Function of Current');
xlim([min(currArr./1000)*0.95,max(currArr./1000)*1.05]);
ylim([min(normErr)*0.95,max(normErr)*1.05]);
xlabel('Current (kA)');
ylabel('Error');
hold off