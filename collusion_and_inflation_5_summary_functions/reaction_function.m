function reaction_function = reaction_function(z1_state,z2_state,a,b,c,r,K)
% Reaction function
% Inputs: 
%   z2 the real prices in logs of the competitor
%   z1 my real price. I need this to determine the current amount of
%   profits, after that z1 is not used
%   (z1, z2) are the current states 
%   a, b, c parameters of the profit function 
%   r, K interest rate and menu cost
% Output: 
%   S(z2), my reaction to z2 

% Left hand side 
profit = quadratic_profit(z1_state,z2_state,a,b,c) / 2;
% C = profit - r*K;

% Solution 
% reaction_function = fsolve (@(x) profit - quadratic_profit(x,z2_state,a,b,c) / 2 + r*K, 1,...
    % optimset('Display', 'off'));

% Parameters of the quadratic 
aa = b; 
bb = -(a+c*z2_state); 
cc = 2*profit - 2*r*K+b*z2_state^2-a*z2_state; 

reaction_function = (-bb + sqrt(bb^2 - 4 * aa * cc)) / (2*aa); 

%reaction_function = (-(4+z2) - sqrt((4+z2).^2 + 4 .* 2 .* (4*z2 - 2*z2.^2 - 8*C))) ...
    %/ (2 * (-2)); 

end