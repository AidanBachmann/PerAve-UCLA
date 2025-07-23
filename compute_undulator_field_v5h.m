% Compute_undulator_field
% Works with perave_core_v5

startindex=max(1,floor(param.z0/param.stepsize));

if param.tapering == 0
    res_phase = zeros(1,param.Nsnap);
elseif param.tapering == 2
    zpos = [1:param.Nsnap]*param.stepsize;
    res_phase = zeros(1,param.Nsnap);
    for ij = 1:param.Nsnap
        if zpos(ij) > param.z0 % Suppress tapering for z < z0
            res_phase(ij) = 1;
        end
    end
else
    res_phase(1:startindex) = 0;
    res_phase(startindex:param.Nsnap) = param.psir;
end
Kz(1) = param.K;