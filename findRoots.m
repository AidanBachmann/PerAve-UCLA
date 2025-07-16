%% Solve cubic equation to find theoretical value of exponential gain length
function [lambda] = findRoots(rho,sigma,delta)
    a = 1; % Coefficients of the cubic polynomial of the form a*lambda^3 + b*lambda^2 + c*lambda + d
    b = (sigma.^2)/2 - delta;
    c = 2*rho - (sigma.^2)/rho;
    d = delta*((sigma.^2)/rho) - (sigma.^4)/(2*rho) + rho*(sigma.^2) + 1;
    rts = imag(roots([a b c d])); % Compute roots, take imaginary part
    lambda = abs(min(rts)); % Find negative root corresponding to exponential growth, take the norm
end