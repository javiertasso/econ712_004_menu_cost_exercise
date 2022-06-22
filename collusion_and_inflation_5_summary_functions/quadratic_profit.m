function quadratic_profit = quadratic_profit(z1,z2,a,b,c)
% Quadratic profit for the multiproduct S&W
% Inputs: 
%   z1 and z2 the real prices in logs
%   a, b, c parameters of the profit function 


quadratic_profit = a * (z1+z2) - b * (z1.^2+z2.^2) + c * z1 .* z2; 

end