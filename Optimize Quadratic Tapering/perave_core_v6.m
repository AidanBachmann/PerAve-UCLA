function [power,radfield,gammap,thetap,profile_l,profile_b,bunch,gammares] = perave_core_v6(param,e0,c,firstpass,Kz,res_phase)
    % This version is for a KMR style 1-D TESSA where you choose the undulator
    % K iteratively by specifying the resonant phase
    %% initialize phase space (Quiet - start problem )
    t_init= tic; % Initial system time
    Np = param.Np; % Grab number of particles
    nbins = 32; % Binning for particles?
    mpart = Np/nbins;
    n_electron = param.I*param.lambda0*param.zsep/e0/c;
    p1 = zeros(Np,1);
    
    tslice = (1:param.nslices)*param.lambda0*param.zsep/c;
    
    if (param.beamdistribution == 1)
        profile_b = exp(-(tslice-tslice(end)/2).^2/2/param.sigma_t^2);
    else
        profile_b(1:param.nslices) = 0;
        profile_b(abs(tslice-tslice(end)/2)<param.sigma_t) = 1;
    end
    
    if (param.laserdistribution == 1)
        profile_l = exp(-(tslice-param.slippage).^2/2/param.sigma_l^2);
    else
        profile_l(1:param.nslices) = 0;
        profile_l(abs(tslice-param.slippage)<param.sigma_l) = 1;
    end
        
    %profile_b = heaviside(tslice-tslice(end)/2+param.sigma_t).*(1-heaviside(tslice-tslice(end)/2-param.sigma_t));
    
    radfield=ones(param.Nsnap,param.nslices)*param.E0;
    radfield(1,:) = profile_l*param.E0;
    
    if ~firstpass
        radfield(1,:) = oldfield;
    end
    
    thetap = zeros(param.Nsnap,param.nslices,Np);
    gammap=zeros(param.Nsnap,param.nslices,Np);
    
    for islice = 1:param.nslices
    X0 = hammersley(2,Np);
    gammap(1,islice,:) = param.gamma0+param.deltagamma*X0(1,:);
    auxtheta1 = hammersley(1,mpart)'*2*pi/nbins-pi;
    for jbin = 1:nbins
        for ipart = 1:mpart
            thetap(1,islice,ipart+(jbin-1)*mpart)=auxtheta1(ipart)+2*(jbin-1)*pi/nbins;
        end
    end
    
    if(param.shotnoise)
        an = 2*sqrt(-log(rand(1))/n_electron);    
        phin = rand(1)*2*pi;
        for ipart = 1:Np;
        thetap(1,islice,ipart) = thetap(1,islice,ipart)-an*sin(thetap(1,islice,ipart)+phin);
        end    
    end
    
    if (param.prebunching ==1 )
        thetap(1,islice,:) = thetap(1,islice,:)-2.*param.bunch*sin(thetap(1,islice,:)+param.bunchphase);
    end
    if (param.prebunching < 0)
           thetab  = squeeze(thetap(1,islice,:));
            gammab = squeeze(gammap(1,islice,:));
            [thetab,gammab] = buncher(thetab,gammab,param.buncherAmp);
            for ipart = 1:Np;
            thetap(1,islice,ipart) = thetab(ipart) + param.bunchphase;
            gammap(1,islice,ipart) = gammab(ipart);
            end
    end
    
    bunching(islice) = (sum(exp(1i.*thetap(1,islice,:))/Np));
    end
    disp(['Initialization time = ',num2str(toc(t_init)),' sec'])
    
    %% Solve system of equations
    
    % define initial conditions
    %options = odeset('RelTol', param.accuracy,'OutputFcn',@odeplot,'OutputSel',[1]); %,'Stats','on');   % set solver options
    %options = odeset('RelTol', param.accuracy); %,'Stats','on'); 
    res_step = param.und_periods*param.lambdau/param.Nsnap;   
    
    total_simtime = 0;
    hl = 0;
    z(1) = 0;
    gammares(1) = sqrt (param.lambdau.*(1+Kz(1)^2) / 2/param.lambda0);            
    param.stepsize = param.lambdau*param.delz;
    % Constant for the resonant phase based tapering   
    const_resp=1/param.chi2*(param.lambdau/2/param.lambda0);
    slip = 0;
    if(param.itdp)
    for ij = 1:param.Nsnap-1  % takes Nsnap snapshots along length of undulator
        tstart = tic;
    % firstslice=ij;  
        for islice = 1:param.nslices       
         gammaf = squeeze(gammap(ij,islice,:));
         thetaf = squeeze(thetap(ij,islice,:));
         E_q0 = radfield(ij,islice);
         param.chi1 = mu0*c/2*param.I*profile_b(islice)/param.A_e;
    % RK4th order integration     
         phasespaceold = [thetaf,gammaf];
         evaluesold = E_q0;
    
         [phasespacenew,evaluesnew]=push_FEL_particles_RK4(phasespaceold,evaluesold,param,Kz(ij));       
    
         thetap(ij+1,islice,:) = phasespacenew(:,1);
         gammap(ij+1,islice,:) = phasespacenew(:,2);
         radfield(ij+1,islice) = evaluesnew;
        end
    % Slippage of the radiation field          
                %Fundamental
                if(mod(ij,param.zsep) == param.zsep-1)
                B = radfield(ij+1,:).';
                B = circshift(B,1);
                radfield(ij+1,:) = B;
    
                if (~firstpass)
                    radfield(ij+1,1)= 0;
                else
                    radfield(ij+1,1) = param.E0*profile_l(1);
                end
                
                end
                
                % Compute bunching
                bunch(ij,:)=(mean(exp(1j.*thetap(ij,:,:)),3)); % shape(thetap) = (Nsnap,nslices,Np)
                
                % Compute undulator field at next step (constant res phase)
                if (firstpass)
                     switch tapering_strength
                         case 0 
                             Klz = max(abs(radfield(1,:)));
                         case 1    
                             Klz = max(abs(radfield(ij,:)));
                         case 2
                             Klz = mean(abs(radfield(ij,1:param.nslices)),2);
                     end
    %         Klz = sum(abs(radfield(ij,:)).*profile_b(:)')/sum(profile_b)/8
                
                Kz(ij+1)=Kz(ij)-param.stepsize/const_resp*Klz.*sin(res_phase(ij));
                end
                
                gammares(ij+1) = sqrt (param.lambdau.*(1+Kz(ij).^2) / 2/param.lambda0);            
                telapsed = toc(tstart);
                total_simtime = total_simtime+telapsed;
        
        formatSpec = '%.3f sec from z = %.3f to z = %.3f, total length = %.3f \n';
        fprintf(formatSpec, telapsed, res_step*(ij-1), res_step*ij, param.und_periods*param.lambdau);
    end
    else
        zpos = [1:param.Nsnap]*param.stepsize;
        %% Time independent
        deltagammamax(1) = 1;
        for ij = 1:param.Nsnap-1  % takes Nsnap snapshots along length of undulator
         gammaf = squeeze(gammap(ij,1,:));
         thetaf = squeeze(thetap(ij,1,:));
         E_q0 = radfield(ij,1);   
         
         % RK4th order integration     
         phasespaceold = [thetaf,gammaf];
         evaluesold = E_q0;
         [phasespacenew,evaluesnew]=push_FEL_particles_RK4(phasespaceold,evaluesold,param,Kz(ij));       
         thetap(ij+1,1,:) = phasespacenew(:,1);
         gammap(ij+1,1,:) = phasespacenew(:,2);
         radfield(ij+1,1) = evaluesnew;          
        
         % Compute undulator field at next step (constant res phase)
                Kz(ij+1) = Kz(ij)-param.stepsize/const_resp*abs(radfield(ij,:)).*sin(res_phase(ij));
                %Kz(ij+1) = param.K*(1 + (param.a0*(zpos(ij) - param.zs) + param.b0*(zpos(ij) - param.zs).^2)*res_phase(ij)); % Quadratic tapering
                gammares(ij+1) = sqrt (param.lambdau.*(1+Kz(ij).^2) / 2/param.lambda0);
                bunch(ij)=(mean(exp(1j.*thetap(ij,1,:)),3));
                param.chi = (Kz(ij+1).^2)/(2*(1 + Kz(ij+1).^2)); % Update chi parameter
          % Particle detrap when deltagamma decreases
                if (ij>40000 & param.changeresphase)
                    bukh = @(phi) sqrt(cos(phi)-(pi/2-phi)*sin(phi));
                    deltagamma = sqrt(Kz(ij+1).*mean(abs(radfield(ij,:)),2))*bukh(res_phase(ij));
                    deltagammamax = max(deltagammamax,deltagamma);
                     if ( deltagamma < deltagammamax)
                         newphase = fsolve(@(phi) sqrt(Kz(ij+1).*mean(abs(radfield(ij,:)),2))*bukh(phi)-deltagammamax, res_phase(ij));  
    %                    res_phase(ij) = newphase-0.2;
    %                    Kz(ij+1)=Kz(ij)-param.stepsize/const_resp*mean(abs(radfield(ij,:)),2).*sin(res_phase(ij));
                     end
     %           res_phase(ij+1) = pi/2 - 2*atan( (Kz( rpindex)*mean(abs(radfield(rpindex,:)),2) ./Kz(ij)./mean(abs(radfield(ij,:)),2))^0.25.*tan(pi/4-param.psir/2));
                end
                
        end
    end
    
    %% Remove slices within one total slippage length
    if(param.itdp)
    radfield(:,1:Nslip)=[];
    gammap(:,1:Nslip,:)=[];
    thetap(:,1:Nslip,:)=[];
    profile_l(1:Nslip)=[];
    profile_b(1:Nslip)=[];
    bunch(:,1:Nslip)=[];
    end
    % Calculate radiation power 
    power(:,:) = abs(radfield(:,:)).^2/377*param.A_e;

end