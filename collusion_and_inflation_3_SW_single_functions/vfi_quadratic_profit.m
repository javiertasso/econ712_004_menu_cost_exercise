function [s,S, z_grid, price_chosen, change_price, V_new] = ...
    vfi_quadratic_profit(g, n_grid, tol_level, ...
    max_n_it, a, b, bbeta, r)
% Performs value function iteration for different inflation rates
%   Inputs
%       n_grid 
%       tol_level
%       max_n_iter  
%       g, inflation rate 
%       (a,b) parameters of the quadratic profit function 
%       (bbeta, r) menu cost, interest rate 

% Create a z grid 
z_grid = transpose(linspace(0,2, n_grid)); 
V_initial = 0 * eye(n_grid,1);  
V_new = 1 * eye(n_grid,1); 
change_price = 0 * eye(n_grid,1); 
price_chosen = 0 * eye(n_grid,1); 
n_it = 0;
diff = 1; 
ddelta = 1 / (1 + r); 

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


end

