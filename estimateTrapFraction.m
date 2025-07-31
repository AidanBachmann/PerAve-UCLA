%% Estimate fraction of total particles trapped in ponderomotive bucket

%{
% Estimate trapped fraction across all snapshots
function ft = estimateTrapFraction(param,thetap)
    ub = 2*pi; % Phase upper bound
    lb = -2*pi; % Phase lower bound
    ft = zeros(param.Nsnap,1); % Trap fraction along undulator
    for ii = 1:param.Nsnap
        counter = 0;
        for jj = 1:param.Np
            if (thetap(ii,1,jj) > lb) && (thetap(ii,1,jj) < ub) % Check if particle is in range of phase values
                counter = counter + 1;
            end
        end
        ft(ii,1) = counter/param.Np; % Compute trapped fraction at iith snapshot
    end
end
%}

% Estimate trapped fraction at single snapshot
function ft = estimateTrapFraction(param,thetap,ii)
    ub = 2*pi; % Phase upper bound
    lb = -2*pi; % Phase lower bound
    counter = 0;
    for jj = 1:param.Np
        if (thetap(ii,1,jj) > lb) && (thetap(ii,1,jj) < ub) % Check if particle is in range of phase values
            counter = counter + 1;
        end
    end
    ft = counter/param.Np; % Compute trapped fraction at iith snapshot
end