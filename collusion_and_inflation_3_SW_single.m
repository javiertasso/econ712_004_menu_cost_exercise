% -------------------------------------------------------------------------
% Collusion and Inflation 3 - S&W Single Product 
% 2022
% Javiera Garcia & Javier Tasso
% Replication of Sheshinski & Weiss: The single product monopoly case
% Try to figure out the value function using a simple value function
% iteration
% Using a quadratic profit function
% -------------------------------------------------------------------------

clearvars
clc 
% cd 'C:\Users\Javier\Documents\Penn\Econ712-Macro-Financial-Markets\Collusion_and_Inflation'
cd;  
addpath([cd '\collusion_and_inflation_3_SW_single_functions'])

% -------------------------------------------------------------------------
% Set values of parameters 

    % Inflation rate 
    g = 0.02;

    % Interest rate 
    r= 0.03; 

    % Discount factor
    ddelta = 1 / (1+r); 
    
    % (real) Menu cost
    bbeta = 0.05;

    % Parameters of profit function
    a = 1; 
    b = 1/2; 

%%
% -------------------------------------------------------------------------

% Create a z grid 
n_grid = 1001;
tol_level = 10^(-5); 
max_n_it = 100000; 
z_grid = transpose(linspace(0,2, n_grid)); 
V_initial = 0 * eye(n_grid,1);  
V_new = 1 * eye(n_grid,1); 
change_price = 0 * eye(n_grid,1); 
price_chosen = 0 * eye(n_grid,1); 
n_it = 0;
diff = 1; 

while diff > tol_level && n_it < max_n_it

    % Construct V_eroded by the inflation rate 
    V_eroded = 0 * eye(n_grid,1);
    for ii = 1:n_grid 
        
        new_price = z_grid(ii,1) - g;
        [min_value, index] = min(abs(new_price - z_grid)); 
    
        if min_value >= g
    
            V_eroded(ii,1) = -10^5;
    
        else
    
            V_eroded(ii,1) = V_initial(index,1);
    
        end 
    
    end
    
    % For loop, for each z on the grid
    for ii = 1:n_grid 
    
        % Doing nothing 
            % Careful: next period I need value function evaluated at z-g
    
            % New price 
            new_price = z_grid(ii,1) - g; 
            [~,new_price_index] = min(abs(z_grid - new_price)); 
    
        V_doing_nothing = a * z_grid(ii,1) - b * z_grid(ii,1)^2 ...
            + ddelta * V_initial(new_price_index,1); 
    
        % Paying fixed cost and changing the price
            % Here I can do a grid search 
            % Careful: next period I need value function evaluated at z-g
    
             
        V_changing_price = max(a * z_grid - b * z_grid.^2 + ddelta * V_eroded) ...
            - bbeta; 
    
        V_new(ii,1) = max(V_doing_nothing, V_changing_price); 

        if V_changing_price > V_doing_nothing

            change_price(ii,1) = 1; 
            [~, index_price_chosen] = max(a * z_grid - b * z_grid.^2 + ddelta * V_eroded); 
            price_chosen(ii,1) = z_grid(index_price_chosen,1); 

        else

            change_price(ii,1) = 0; 
            price_chosen(ii,1) = z_grid(ii,1); 

        end
    
    end 

    diff = max(abs(V_new - V_initial));
    n_it = n_it + 1; 
    V_initial = V_new; 

end

change_price_temp = change_price(2:end,1); 
[~, index_min] = min(change_price(1:(end-1),1) - change_price_temp); 
S = z_grid(index_min); 
[~, index_max] = max(change_price(1:(end-1),1) - change_price_temp); 
s = z_grid(index_max); 
clear change_price_temp intex_min index_max 
display(s)
display(S)

figure(1)
plot(z_grid(floor(n_grid/4):floor(3*n_grid/4)), ...
    (1-ddelta) * V_new(floor(n_grid/4):floor(3*n_grid/4)), ...
    'LineWidth',2) 
hold on
xlabel('z')
ylabel('V(z) - Value Function')
title('Value function of single firm - single product case')
saveas(gcf, 'collusion_and_inflation_3_value_fun.png')
close(figure(1))

figure(2)
plot(z_grid(floor(n_grid/4):floor(3*n_grid/4)), ...
    change_price(floor(n_grid/4):floor(3*n_grid/4)), ...
    'LineWidth',2) 
hold on
xlabel('z')
ylabel('p - Policy Function')
title('Policy function of single firm - single product case')
saveas(gcf, 'collusion_and_inflation_3_policy_fun.png')
close(figure(2))

figure(3)
plot(z_grid(floor(n_grid/4):floor(3*n_grid/4)), ...
    price_chosen(floor(n_grid/4):floor(3*n_grid/4)), ...
    'LineWidth',2) 
hold on
xlabel('z')
ylabel('Price - Optimal Choice')
title('Price chosen of single firm - single product case')
saveas(gcf, 'collusion_and_inflation_3_price_chosen.png')
close(figure(3))

%%
% -------------------------------------------------------------------------
clearvars 

% Parameters
r = 0.03;
a = 1;
b = 0.5;
bbeta = 0.05; 

% For different inflation rates 
N = 11; 
g_grid = transpose(linspace(0.005,0.1,N));
s = 0 * eye(N,1);
S = 0 * eye(N,1); 

for ii = 1:N

    [s(ii,1),S(ii,1),~,~,~] = vfi_quadratic_profit(g_grid(ii,1),1001,10^(-5), ...
        100000, a, b, bbeta, r);

end

figure(4)
plot(g_grid, s, 'LineWidth',2)
hold on
plot(g_grid, S, 'LineWidth',2)
xlabel('Inflation')
ylabel('(s,S)')
title('Threshold for different inflation rates')
saveas(gcf, 'collusion_and_inflation_3_sS_inflation.png')
close(figure(4))