% Compute_undulator_field
% Works with perave_core_v5

startindex=max(1,floor(param.z0/param.stepsize));

if param.tapering==0
    res_phase=zeros(1,param.Nsnap);
else
    res_phase(1:startindex)=0;
    res_phase(startindex:param.Nsnap)=param.psir;
end
Kz(1)=param.K;
