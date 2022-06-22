function [real_price] = real_price(g,q,t)
% (log of) Real price
%   Inputs: 
%       g, inflation rate
%       q, nominal price 
%       t, time
%   Output: 
%       w_t, log of real price 

real_price = log(q) - g * t; 

end