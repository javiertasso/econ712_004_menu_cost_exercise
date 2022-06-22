function [V_tilde] = V_tilde(a, b, g, bbeta, r)
% V_tilde computes discounted profits between 0 and T, with corresponding
% lower bound a and upper bound b. Important: V_tilde assumes the other
% firm charges the same price, so profits are split evenly
%   Inputs: 
%       a, lower bound, corresponding to t=0
%       b, upper bound, corresponding to t=T
%       g, inflation rate
%       bbeta, menu cost 
%       r, discount rate
%   Output: 
%       V_tilde, discounted profit


fun = @(x) ((exp(x) .* (1-exp(x)))/2) .* exp(-(r/g)*(a-x)); 
I = integral(fun,a,b); 

V_tilde = (1/g) * I - bbeta * exp(-(r/g)*(a-b)); 

end