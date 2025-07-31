function [phaseb,energyb] = buncher(thetab,gammab,amp)
amp1 = amp/4.7;
BR56 = pi /2/(amp+3);
BR561 = pi/(amp1+3);

figure(3)
gamma0 = mean(gammab);
gammaspread = std(gammab);
plot(thetab,gammab)
gammarel = (gammab - gamma0) / gammaspread;
gammarel = gammarel - amp1*sin(thetab);
figure(6)
plot(thetab,gammarel)
phaseb = thetab+gammarel*BR561;
figure(4)
plot(phaseb,gammarel)
gammarel = gammarel - amp*sin(phaseb);
phaseb = phaseb+gammarel*BR56;
energyb = gammarel*gammaspread+gamma0;
figure(5)
plot(phaseb,energyb)
end
