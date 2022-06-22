function [nominal_price] = nominal_price(g,w,t)
% Nominal price
%   Inputs: 
%       g, inflation rate
%       w, (log of) real price 
%       t, time
%   Output: 
%       q_t, nominal price

nominal_price = exp(w+g*t); 

end