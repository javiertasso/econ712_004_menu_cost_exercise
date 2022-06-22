% -------------------------------------------------------------------------
% Collusion and Inflation 2 - S&W Multiproduct 
% 2022
% Javiera Garcia & Javier Tasso
% Replication of Sheshinski & Weiss: The multiproduct monopoly case
% Using a quadratic profit function
% -------------------------------------------------------------------------

clearvars
clc 
% cd 'C:\Users\Javier\Documents\Penn\Econ712-Macro-Financial-Markets\Collusion_and_Inflation'
cd
addpath([cd '\collusion_and_inflation_2_SW_multi_functions'])

%--------------------------------------------------------------------------
% Set values of parameters 

    % Inflation rate 
    g = 0.02;

    % Interest rate 
    r= 0.03; 
    
    % (real) Menu cost
    bbeta = 0.05;

    % Parameters of profit function
    a = 1; 
    c = 1/2; 
    b = 1/2; 

% Profit function is given by the following expression
    % F(z1,z2) = a(z1+z2) - b(z1^2+z2^2) + cz1z2
        % b > 0
        % 4b^2-c^2 > 0
        % c > 0, complementarity between real prices (in logs) 
        % z are real prices in logs
    % Try the function 
    display(quadratic_profit(1,1,a,b,c))
    display(quadratic_profit(2,2,a,b,c))
    display(quadratic_profit(2.5,2.5,a,b,c))

% Try the sss function
display(sss(4,1,r,g,bbeta,a,b,c))
display(sss(5,1,r,g,bbeta,a,b,c))
display(sss(40,1,r,g,bbeta,a,b,c))

% Find (S,epsilon) for this case 
[temp, ~] = fsolve(@(x) sss(x(1),x(2),r,g,bbeta,a,b,c),[1,1], optimset('Display','off'));
S_example = temp(1,1);
epsilon_example = temp(1,2); 
clear temp 

% Figure out little s now 
s_example = S_example / exp(g*epsilon_example);
display(S_example)
display(s_example)
display(epsilon_example)
clear s_example S_example epsilon_example

% How do thresholds behave for different inflation rates 
n_grid = 101; 
inflation_grid = transpose(linspace(10^(-4),0.1,n_grid)); 
S = eye(n_grid, 1);
epsilon = eye(n_grid, 1); 

for ii = 1:n_grid

[temp, ~] = fsolve(@(x) sss(x(1),x(2),r,inflation_grid(ii,1),bbeta,a,b,c),...
    [1,1], optimset('Display','off'));
S(ii,1) = temp(1,1); 
epsilon(ii,1) = temp(1,2);

end

s = S ./ exp(inflation_grid .* epsilon);

clear temp ii 

figure(1)
plot(inflation_grid, s, 'linewidth',2)
hold on
plot(inflation_grid, S, 'linewidth',2)
xlabel('Inflation')
ylabel('(s,S) thresholds')
title('(s,S) for different inflation rates')
subtitle('Syncronized Steady State') 
saveas(gcf, 'collusion_and_inflation_2_sS_bands_synchronized_ss.png')
close(figure(1))

figure(2)
plot(inflation_grid, s, 'linewidth',2)
hold on
xlabel('Inflation')
ylabel('Time without changing prices')
title('Time without changing the price for different inflation rates')
subtitle('Syncronized Steady State') 
saveas(gcf, 'collusion_and_inflation_2_epsilon_synchronized_ss.png')
close(figure(2))




