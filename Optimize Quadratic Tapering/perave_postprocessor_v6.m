%% Perave_postprocessor
close all
kw=2*pi/param.lambdau;
hbar=6.582e-16;
zpos= [1:param.Nsnap]*param.stepsize;

%% Spectrum as a function of z
fundpower=[];
sidebandpower=[];
if param.itdp
omegamin=-10e-4; omegamax=10e-4;
h=figure(4);
set(h, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
for n=1:param.Nsnap
[powerspec,omega]=spectrum_calc((radfield(n,:)),param.lambda0,param.zsep);
sidebandindex=omega>omegamin & omega<omegamax;
fundspectrum=powerspec(sidebandindex);
fundpower(n)=trapz(fundspectrum)/trapz(powerspec);
end
figure(4)
subplot(1,2,1)
semilogy(omega/rho1D,abs(powerspec))
xlabel('\delta\omega/\rho\omega ','FontSize',16)
    ylabel('P (\omega) [arb. units]','FontSize',16)    
    xlim([-20,20])    
    set(gca,'FontSize',16)
    legend(sprintf(['z / L_u =',num2str(zpos(n)/lwig)]));
    
subplot(1,2,2)
plot([1:1:size(power,2)]*param.zsep*param.lambda0*1e15/3e8,power(n,:))
xlim([1,size(power,2)]*param.zsep*param.lambda0*1e15/3e8)
xlabel('t [fs]','FontSize',16)
ylabel('Output Radiation Power [W]','FontSize',16)
set(gca,'FontSize',16)
end

%% Radiation Power and spectrum at exit
avgPower = mean(power,2); % Compute average power
if param.P0 < param.Psat
    pks = findpeaks(avgPower);
    if isempty(pks) == true
        idx = find(avgPower==max(avgPower));
    else
        idx = find(avgPower==pks(1)); % Find index where power is maximal
    end
    idxPad = round(param.Nsnap*0.125); % Curve is not linear near max, use to truncate data for fit before max
    upperIdx = idx - idxPad; % Compute upper index for fit
    lowerIdx = round(idxPad*1.5); % Compute lower index for fit
    coeff = polyfit(zpos(lowerIdx:upperIdx),log(avgPower(lowerIdx:upperIdx)),1); % Compute coefficients for linear fit
    gainlen = 1/coeff(1);
    param.gainlen = gainlen; % Store numerical gain length
    powerFit = exp(polyval(coeff,zpos)); % Compute fit
end
if param.suppress_plots == 0
    figure(2)
    title('Simulation Output')
    subplot(2,2,1)
    % Plots
    semilogy(zpos,avgPower,'b');
    hold on
    semilogy(zpos,max(power'),color='r');
    if param.P0 < param.Psat
        semilogy(zpos,powerFit,color='g',LineStyle='--'); % Plot fit
        xline([zpos(lowerIdx) zpos(upperIdx)],'--',color='r');
        legend('Avg','Max',strcat('Fit, Gain Length =  ',num2str(gainlen)),'Fitting Region');
        fprintf('\nTheoretical Gain Length: %f\nNumerically Computed Gain Length: %f\nNormalized Error: %f\n',param.Lgain,gainlen,abs(gainlen - param.Lgain)/param.Lgain);
    else
        legend('Avg','Max');
    end
    xlim([0,zpos(end)]);
    title('Radiation Power along the beam');
end
if param.itdp
subplot(2,2,2)
plot([1:1:size(power,2)]*param.zsep*param.lambda0*1e15/3e8,power(end,:))
hold on
plot([1:1:size(power,2)]*param.zsep*param.lambda0*1e15/3e8,max(power(end,:))*profile_l(:),'-')
plot([1:1:size(power,2)]*param.zsep*param.lambda0*1e15/3e8,max(power(end,:))/2*profile_b(:),'.')
xlim([1,size(power,2)]*param.zsep*param.lambda0*1e15/3e8)
xlabel('t [fs]')
ylabel('Power [W]')
legend('Final',sprintf(['Norm Initial ',num2str(P0/1e9),' GW']), sprintf(['Current profile ',num2str(param.I/1e3,3),' kA']),'Location','southeast')

[powerspec,omega]=spectrum_calc(radfield(end,:),param.lambda0,param.zsep);
subplot(2,2,3)
semilogy((omega+1)*hbar*2*pi*c/param.lambda0,powerspec,'b');    
xlim([omega(1)+1,omega(end)+1]*hbar*2*pi*c/param.lambda0)    
    xlabel('Photon Energy [eV]')
    ylabel('P (\omega) [arb. units]')
    title('Output Spectrum')
end

%% Bunching factor and energy loss 
%{
figure(2)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);

subplot(2,3,4)
plot([1:1:param.Nsnap-1]*param.stepsize,mean(abs(bunch),2))
xlim([0,param.Nsnap*param.stepsize])
xlabel('z')
ylabel('Bunching Factor')

subplot(2,3,5)
for ij = 1:param.Nsnap-1
    bucket_amp_rel =    2*sqrt(Kz(ij)*param.chi2*param.E0*param.lambda0/2/pi/(1+Kz(ij)^2));
    idx = squeeze(gammares(ij)*(1+4*bucket_height*bucket_amp_rel)>gammap(ij,:,:));
    trap(ij) = sum(sum(idx,2)/param.Np.*profile_b')/sum(profile_b);
end
xlabel('z')
ylabel('trapped fraction')
plot([1:1:param.Nsnap-1]*param.stepsize,trap)
%}

for ij=1:param.Nsnap
meanenergy(ij)=sum(mean(gammap(ij,:,:),3).*profile_b)/sum(profile_b);
end

% 
if param.suppress_plots == 0
    subplot(2,2,4)
    plot([1:1:param.Nsnap]*param.stepsize,meanenergy)
    xlabel('z')
    ylabel('\gamma')
    xlim([0,param.Nsnap*param.stepsize])
end

if param.itdp
end
%% Sven's plots
if param.itdp
    figure(3)
    subplot(1,3,1)
    contour(abs(bunch))
    subplot(1,3,2)
    contour(angle(bunch))
    subplot(1,3,3)
    contour(power)  
end

%% Energy calculations 
if(param.itdp)
    pulselength = fwhm([1:size(power,2)]*param.zsep*param.lambda0*1e15/3e8,smooth(power(end,:),15))
else
    pulselength = param.sigma_t*2.35;
end

Ebeam = meanenergy(1)*param.I*sum(profile_b)*param.lambda0*param.zsep/c*511000;
Erad = (sum(power(end,:))-sum(power(1,:)))*param.lambda0*param.zsep/c;
Eloss = (meanenergy(end)-meanenergy(1))*param.I*sum(profile_b)*param.lambda0*param.zsep/c*511000;
Efficiency = Erad/Ebeam
Econservation = ((meanenergy(end)-meanenergy(1))*param.I*sum(profile_b)*param.lambda0*param.zsep/c*511000+Erad)/Ebeam
%% % Phasespace movie
if param.phasespacemovie
    filename='particle_movie.gif';
    figure(10)
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    for i=1:param.Nsnap
        fieldphase(i)=mean(angle(radfield(i,:)));
        
        tp=squeeze(thetap(i,:,:))+fieldphase(i)+pi/2;% You can add the phase of the field if you want
        gp=squeeze(gammap(i,:,:));    
        
        tresh=reshape(tp,[1,size(tp,1)*size(tp,2)]);
        gresh=reshape(gp,[1,size(gp,1)*size(gp,2)]);
        
        if param.itdp
            plot((mod(tresh,2*pi)-pi)./pi,(gresh./(meanenergy(1))-1),'*k','MarkerSize',1)
        else     
           plot((mod(tp+2*pi,2*pi)-pi)./pi,(gp./(meanenergy(1))-1),'.k','MarkerSize',3)
        end
        set(gca,'FontSize',20)
        xlabel('\Psi/pi');ylabel('\Delta \gamma/\gamma_0');
        pause(0.1)
    end
end
    
