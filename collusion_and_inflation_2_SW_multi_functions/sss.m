function sss = sss(S,epsilon,r,g,bbeta,a,b,c)
% sss for syncronized steady state
%   Detailed explanation goes here

% Function to integrate 
fun = @(x) exp(-r*x) .* (a*(2*S-2*g*x) - b*(2*(S-g*x).^2) + c*(S-g*x).^2); 

% Integral 
II = integral(fun, 0, epsilon);

eq1 = quadratic_profit(S,S,a,b,c) - r / (1-exp(-r*epsilon)) * II + 2 * bbeta * exp(-r*epsilon);
eq2 = quadratic_profit(S-g*epsilon,S-g*epsilon,a,b,c) - quadratic_profit(S,S,a,b,c) ...
    + 2*r*bbeta;
sss(1,1) = eq1;
sss(2,1) = eq2;
end