function [max_V, z1_star, z2_star, z1_low, z1_high, z2_low, z2_high, V_fun, ...
    z1_grid, z2_grid, V_opt] = ...
    VFI_competitive_case_1(n_grid, tol_level, max_n_it, a, b, c, ...
    r, K, g, V_initial_guess)
% Performs VFI for the competitive case when q=1
%   Inputs
%       n_grid, tol_level, max_n_it to initialize the algorithm
%       a, b, c for the profit function
%       r, K and g to determine beta, menu cost and inflation rate
%   Outputs
%       max_V, z1_star z2_star
%       thresholds of inaction
%       the whole value function in a matrix + grids
%       V_opt

% VFI Competitive case with q = 1 -----------------------------------------

% Generate a grid
%n_grid = 101;
%tol_level = 10^(-2);
%max_n_it = 1000; 
z1_grid = transpose(linspace(1,2,n_grid)); 
z2_grid = transpose(linspace(1,2,n_grid)); 

if isempty(V_initial_guess)

    V_initial = 15 * ones(n_grid, n_grid); 

else

    V_initial = V_initial_guess;

end

V_new = 0 * eye(n_grid,n_grid); 
n_it = 0;
diff = 1; 
bbeta = 1 / (1 + r);
change_price = 0 * eye(n_grid, n_grid); 
price_1_chosen = 0 * eye(n_grid, n_grid); 
price_2_chosen = 0 * eye(n_grid, n_grid); 

% The following is going to be inside a while 
while diff > tol_level && n_it < max_n_it
    
    % Construct V_eroded by the inflation rate 
    V_eroded = 0 * eye(n_grid, n_grid); 

    for ii = 1:n_grid

       new_price_1 = z1_grid(ii,1) - g; 
       [~, index_min_1] = min(abs(new_price_1 - z1_grid)); 

       for jj = 1:n_grid 

            new_price_2 = reaction_function(z1_grid(ii,1), z2_grid(jj,1), ...
                a, b, c, r, K) - g;
            [~, index_min_2] = min(abs(new_price_2 - z2_grid));
            
            V_eroded(ii,jj) = V_initial(index_min_1, index_min_2); 
        
       end

    end
    
    V_eroded(1,:) = -10^(5);
    V_eroded(:,1) = -10^(5);  

     

    % For loop 

    for ii = 1:n_grid

        for jj = 1:n_grid 

            % Do nothing 
            new_price_1 = z1_grid(ii,1) - g;
            [~,new_price_1_index] = min(abs(z1_grid - new_price_1)); 
            new_price_2 = reaction_function(z1_grid(ii,1), z2_grid(jj,1), ...
                a,b,c,r,K) - g;
            [~,new_price_2_index] = min(abs(z2_grid - new_price_2));
            reaction_price_2 = reaction_function(z1_grid(ii,1), z2_grid(jj,1), ...
                a,b,c,r,K); 
            [~,reaction_price_index] = min(abs(z2_grid - reaction_price_2));
  
            V_doing_nothing = quadratic_profit(z1_grid(ii,1), ...
                z2_grid(reaction_price_index,1), a, b, c) / 2 ...
                + bbeta * V_initial(new_price_1_index, new_price_2_index);
             
            % Pay fixed cost and change prices 
                
                % First be careful with the reaction function
                % Here I create a reaction vector given z2
                reaction_vector_temp = 0 * eye(n_grid, 1);
                reaction_vector_indi = 0 * eye(n_grid, 1);
                reaction_vector_indi_fut = 0 * eye(n_grid, 1); 

                for kk = 1:n_grid 

                    reaction_vector_temp(kk,1) = reaction_function(z1_grid(kk,1), ...
                        z2_grid(jj,1),a,b,c,r,K); 
                    [~, reaction_vector_indi(kk,1)] = min(abs(z2_grid - ...
                        reaction_vector_temp(kk,1))); 
                    [~, reaction_vector_indi_fut(kk,1)] = min(abs(z2_grid - ...
                        reaction_vector_temp(kk,1) + g)); 

                end
                
                % Now I want to find z1 that maximizes my profit 
                % To do this, first try to match the indices 
                V_changing_price_temp = 0 * eye(n_grid, 1); 

                for kk = 1:n_grid 

                    [~,index_fut] = min(abs(z1_grid - z1_grid(kk,1) + g)); 

                    V_changing_price_temp(kk,1) = quadratic_profit(z1_grid(kk,1), ...
                        reaction_vector_temp(kk,1), a, b, c) / 2 + bbeta * ...
                        V_eroded(index_fut,reaction_vector_indi_fut(kk,1)); 

                end

                V_changing_price = max(V_changing_price_temp) - K;               
                [~, p1_chosen_index] = max(V_changing_price_temp);
                p2_chosen_index = reaction_vector_indi(p1_chosen_index); 
 
                if V_changing_price > V_doing_nothing

                    change_price(ii,jj) = 1;
                    price_1_chosen(ii,jj) = z1_grid(p1_chosen_index,1);
                    price_2_chosen(ii,jj) = z2_grid(p2_chosen_index,1); 

                else

                    change_price(ii,jj) = 0;
                    price_1_chosen(ii,jj) = z1_grid(ii,1);
                    price_2_chosen(ii,jj) = reaction_vector_indi(ii,1); 

                end

                V_new(ii,jj) = max(V_doing_nothing, V_changing_price); 
               
        end

    end

     

    diff = max(max(abs(V_new - V_initial)));
    n_it = n_it + 1; 
    V_initial = V_new; 

end

% Outputs

% Max of the value function
max_V = (1-bbeta) * max(max(V_new));

% Where the max is located 
[z1_star_index, z2_star_index] = find(V_new == max_V/(1-bbeta)); 
z1_star = z1_grid(z1_star_index);
z2_star = z2_grid(z2_star_index);

% Thresholds 
V_temp_1 = V_new(:,z2_star_index);
V_temp_1_l = [V_temp_1(1,1); V_temp_1(1:(end-1),1)];
V_temp_1_f = [V_temp_1(2:end, 1); V_temp_1(n_grid,1)]; 

jump = 0; 
ii = 0; 
while jump == 0

    ii = ii + 1; 
    jump = V_temp_1(ii) - V_temp_1_l(ii); 

end
z1_low = z1_grid(ii); 

jump = 0; 
ii = n_grid;
while jump == 0

    ii = ii-1;
    jump = V_temp_1(ii) - V_temp_1_f(ii); 

end
z1_high = z1_grid(ii); 

V_temp_2 = transpose(V_new(z1_star_index,:));
V_temp_2_l = [V_temp_2(1,1); V_temp_2(1:(end-1),1)];
V_temp_2_f = [V_temp_2(2:end, 1); V_temp_2(n_grid,1)]; 

jump = 0; 
ii = 0; 
while jump == 0

    ii = ii + 1; 
    jump = V_temp_2(ii) - V_temp_2_l(ii); 

end
z2_low = z2_grid(ii); 

jump = 0; 
ii = n_grid;
while jump == 0

    ii = ii-1;
    jump = V_temp_2(ii) - V_temp_2_f(ii); 

end
z2_high = z2_grid(ii); 

clear z1_star_index z2_star_index ii jump V_temp_1 V_temp_1_f V_temp_1_l ...
    V_temp_2 V_temp_2_l V_temp_2_f

% Value function
V_fun = (1-bbeta) * V_new; 

% V optimal at the choice z1h, z2h
V_opt = V_fun(z1_grid == z1_high, z2_grid == z2_high); 

end