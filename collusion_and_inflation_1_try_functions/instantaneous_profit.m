function instantaneous_profit = instantaneous_profit(w,w_others)
% Instantaneous profit of firm i, given its price + competitors 
%   Inputs: 
%       Your real price at date t, that is, wt
%       Vector of real prices of competitors 
%       Time 
%   Output: 
%       Instantaneous profit at that moment 
%   Comments: 
%       Let w be the real price in log
%       Then z=exp(w) is the real price
%       This is the profit F(z) = z * (1-z) 
%           This is demand Q=(1-p) and marginal cost equal to zero
%       Additionally I want z between 0 and 1 

w_others_min = min(w_others); 

z = exp(w); 

if w < w_others_min

    instantaneous_profit = z * max((1-z),0);

end

if w > w_others_min 

    instantaneous_profit = 0;

end

if w == w_others_min

    instantaneous_profit = (z * max((1-z),0)) / (length(w_others) + 1); 

end

end