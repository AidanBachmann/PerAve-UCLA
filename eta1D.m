% Compute the theoretical 1D relative energy loss along the undulator.
function [eta_1d] = eta1D(param,bunch,zi)
    JJ = besselj(0,param.chi) - besselj(1,param.chi);
    eta_1d = ((param.chi2/param.gamma0)*((param.K*JJ*bunch)/(2*param.gamma0)))*(param.E0*zi + param.chistar*((param.K*JJ*bunch)/(2*param.gamma0))*zi.^2);
end