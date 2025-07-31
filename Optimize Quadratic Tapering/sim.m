function [pmax] = sim(x0)
    global c mu0 e0 me eps0 IA Z0
    
    c = 2.99792458.*10^8;                                            % speed of light
    e0 = 1.60217657e-19;                                             % electron charge
    me = 9.10938291e-31;                                             % electron mass
    eps0 = 8.85418782e-12;                                           % eps_0
    mu0 = 1.256637e-6;                                               % mu_0
    IA = 17045;                                                      % Alfven current
    Z0 = 376.73;   
        
    param.sigma_t = 0.9e-12;
    param.charge = 2255.965e-12; % Charge
    param.use3Dcorrection  = 1;
    param.beamdistribution = 1;       % Using GENESIS flag: 2-uniform 1-gaussian
    param.laserdistribution = 1;         % Using GENESIS flag: 2-uniform 1-gaussian
    recirculate = 0;
    t1 = tic;

    Perave_User_Input_osc;
    param.a0 = x0(1); % Set parameters for quadratic tapering
    param.b0 = x0(2);
    param.zs = x0(3);

    %% Compute the undulator field
    [Kz,res_phase] = compute_undulator_field_v5h(param);
    
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
    firstpass = 1;
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
    [power,radfield,gammap,thetap,profile_l,profile_b,bunch,gammares] = perave_core_v6(param,e0,c,firstpass,Kz,res_phase);
    disp(['Simulation time = ',num2str(toc(t0)./60),' min'])
    perave_postprocessor_v6   
    rad_vs_und(:,npasses) = sum(power,2)*param.lambda0*param.zsep/c;
    rad_vs_beam(:,npasses) = power(end,:);
    Eff(npasses) = Efficiency;
    PL(npasses) = pulselength;
    oldfield(1:param.nslices) = 0;
    
    if cavitydetuning>0
    oldfield(1,cavitydetuning+1:cavitydetuning+size(radfield,2)) = radfield(end,:)*sqrt(transmission);
    else
    oldfield(1,1:1+cavitydetuning+size(radfield,2)) = radfield(end,-cavitydetuning:end)*sqrt(transmission);    
    end
    pause(0.5)

    %%
    %figure(8)
    %subplot(1,2,1)
    %plot(abs(fftshift(fft(oldfield)).*filter3));
    %subplot(1,2,2)
    filterfield = ifft(ifftshift(fftshift(fft(oldfield) ).*filter3));
    %plot(power(end,:),'k')
    %hold on
    %plot(abs(filterfield).^2/377*param.A_e,'g')
    %plot(abs(oldfield).^2/377*param.A_e,'r')
    %plot(profile_b*max(power(end,:))*0.5,'b')
    
    %hold off
    %pause(0.5)
    oldfield = filterfield;
    firstpass = 0; 
    end
    if param.suppress_plots == 0
        gamma_avg = transpose(mean(gammap,3)); % Average gamma along undulator
        detune = (gamma_avg.^2 - gammares)./(2*(gammares.^2)*param.rho1D); % Detuning along undulator
        figure(5);
        plot([1:1:param.Nsnap]*param.stepsize,detune);
        title('Detuning Along Undulator');
    end
    pmax = -1*max(avgPower);
    pmax
end