function [SW_system] = sheshinsky_weiss(s,S,r,g,bbeta)
% Given s and S, compute the system of equations that I want to solve to
% find S&W thresholds. 
%   Use the following function F(z) = z * (1-z) / 2, that is split monopoly
%   profits by the number of firms + assume linear demanda and zero
%   marginal cost equal between firms. The nice thing about this example is
%   that I can solve the integral by hand 

eq1 = (s * (1-s)) / 2 - (S * (1-S)) / 2 + r*bbeta;
eq2 = (1/2) * ((S^(r/g+1) / (r/g+1) - (2*S^(r/g+2))/(r/g+2)) - ...
    (s^(r/g+1) / (r/g+1) - (2*s^(r/g+2))/(r/g+2))); 

SW_system(1,1) = eq1;
SW_system(2,1) = eq2; 

end