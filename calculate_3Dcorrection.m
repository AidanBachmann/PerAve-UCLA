%% Three Ming Xie parameters
eta_diffr = Lgain / (4*pi*param.sigmax^2/param.lambda0);
eta_em = Lgain / betax *(4*pi*emitx/param.lambda0/param.gamma0); 
eta_es = Lgain *(4*pi/param.lambdau)*param.deltagammarel;

a = [0.45, 0.57, 0.55, 1.6, 3, 2, 0.5, 2.9, 2.4, 51, 0.95, 3, 5.4, 0.7, 1.9, 1140, 2.2, 2.9, 3.2];

eta = a(1)*eta_diffr.^a(2)+a(3)*eta_em^a(4)+a(5)*eta_es^a(6);
eta = eta+a(7)*eta_em^a(8)*eta_es^a(9)+a(10)*eta_diffr^a(11)*eta_es^a(12)+a(13)*eta_diffr^a(14)*eta_em^a(15);
eta = eta+a(16)*eta_diffr^a(17)*eta_em^a(18)*eta_es^a(19);

eta
Lgain3D = Lgain*(1+eta)
