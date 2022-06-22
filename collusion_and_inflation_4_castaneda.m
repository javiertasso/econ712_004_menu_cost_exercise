% -------------------------------------------------------------------------
% Collusion and Inflation 4 - S&W competition case  
% 2022
% Javiera Garcia & Javier Tasso
% Following Castanneda phd thesis (1994)
% Using a quadratic profit function
% -------------------------------------------------------------------------

clearvars
clc 
% cd 'C:\Users\Javier\Documents\Penn\Econ712-Macro-Financial-Markets\Collusion_and_Inflation'
cd
addpath([cd '\collusion_and_inflation_4_castaneda_functions'])

% Parameters
    % (real) Menu cost 
    K = 0.05;

    % Inflation rate
    g = 0.02; 

    % Interest rate 
    r = 0.03; 

    % Parameters of profit function
    a = 1; 
    c = 1/4; 
    b = 1/2; 

% Try profit function 
    display(quadratic_profit(1,1,a,b,c))
    display(quadratic_profit(2,2,a,b,c))
    display(quadratic_profit(2.5,2.5,a,b,c))

% First I want to solve the collusive case 
    % I'm gonna assume I decide to change prices at the same time 

% VFI Colusive case -------------------------------------------------------

% Generate a grid
n_grid = 101;
tol_level = 10^(-2);
max_n_it = 1000; 
z1_grid = transpose(linspace(1,2,n_grid)); 
z2_grid = transpose(linspace(1,2,n_grid)); 
V_initial = 30 * ones(n_grid, n_grid); 
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
       [min_value_1, index_min_1] = min(abs(new_price_1 - z1_grid)); 

       for jj = 1:n_grid 

            new_price_2 = z2_grid(jj,1) - g;
            [min_value_2, index_min_2] = min(abs(new_price_2 - z2_grid));
            
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
            new_price_2 = z2_grid(jj,1) - g;
            [~,new_price_2_index] = min(abs(z2_grid - new_price_2));

            V_doing_nothing = quadratic_profit(z1_grid(ii,1), z2_grid(jj,1), a, b, c)...
                + bbeta * V_initial(new_price_1_index, new_price_2_index);
            
            % Pay fixed cost and change prices 
                % Grid search being careful

                V_temp = 0 * eye(n_grid,1);
                p2_chosen_index_v = 0 * eye(n_grid,1); 

                for kk = 1:n_grid % loop over the first price 

                    [V_temp(kk,1), p2_chosen_index_v(kk,1)] = ...
                        max(quadratic_profit(z1_grid(kk), z2_grid, a, b, c)...
                        + bbeta * transpose(V_eroded(kk,:))); 

                end

                V_changing_price = max(V_temp) - 2*K;
                [~, p1_chosen_index] = max(V_temp);
                p2_chosen_index = p2_chosen_index_v(p1_chosen_index); 

                if V_changing_price > V_doing_nothing

                    change_price(ii,jj) = 1;
                    price_1_chosen(ii,jj) = z1_grid(p1_chosen_index,1);
                    price_2_chosen(ii,jj) = z2_grid(p2_chosen_index,1); 

                else

                    change_price(ii,jj) = 0;
                    price_1_chosen(ii,jj) = z1_grid(ii,1);
                    price_2_chosen(ii,jj) = z2_grid(jj,1); 

                end

                V_new(ii,jj) = max(V_doing_nothing, V_changing_price); 
              
        end

    end

    

    diff = max(max(abs(V_new - V_initial)));
    n_it = n_it + 1; 
    V_initial = V_new; 

end 

% It seems to work 
plot(z2_grid, V_new(floor(n_grid/2),:))
plot(z2_grid, price_2_chosen(floor(n_grid/2),:))
figure(1)
surf(z1_grid, z2_grid, (1-bbeta)*V_new)
saveas(gcf, 'collusion_and_inflation_4_value_fun.png')
close(figure(1))
surf(z1_grid, z2_grid, price_1_chosen)
surf(z1_grid, z2_grid, price_2_chosen)

contour(z1_grid, z2_grid, price_1_chosen,8)

% Plots
figure(1)
plot(z1_grid, (1-bbeta) * V_new(:,floor(n_grid/4)), 'LineWidth',2)
hold on
plot(z1_grid, (1-bbeta) * V_new(:,floor(2*n_grid/4)), 'LineWidth',2)
plot(z1_grid, (1-bbeta) * V_new(:,floor(3*n_grid/4)), 'LineWidth',2)
xlabel('z1')
ylabel('Value Function')
title('Value Function Collusive Case')
legend({strcat('z2= ', num2str(z1_grid(floor(n_grid/4)))), ...
    strcat('z2= ', num2str(z1_grid(floor(n_grid/2)))), ...
    strcat('z2= ', num2str(z1_grid(floor(3*n_grid/4))))}, ...
    'Location', 'northeast')
saveas(gcf, 'collusion_and_inflation_4_value_fun_collusion.png')
close(figure(1))


display((1-bbeta)*max(max(V_new)))

display((1-bbeta)*max(max(V_new)))

% What's next?
    % 1) This is the colusive case. They are going to split this value in
    % two. 
    % 2) I need to figure out how much profits they get in the competitive
    % case
        % This implies I need to understand the reaction function shit 
    % 3) I need to figure out the value of deviating
    % 4) Finally get the discount factor and repeat for another inflation
    % rate
      
%% 
% -------------------------------------------------------------------------
clearvars 
% Parameters
    % (real) Menu cost 
    K = 0.05;

    % Inflation rate
    g = 0.02; 

    % Interest rate 
    r = 0.03; 

    % Parameters of profit function
    a = 1; 
    c = 1/4; 
    b = 1/2; 

% VFI when q=0 ------------------------------------------------------------

% Here I can take z2 as given 
% This seems to work OK

% Generate a grid
n_grid = 101;
tol_level = 10^(-2);
max_n_it = 1000; 
z1_grid = transpose(linspace(1,2,n_grid)); 
z2_grid = transpose(linspace(1,2,n_grid)); 
V_initial = 15 * ones(n_grid, n_grid); 
V_new = 0 * eye(n_grid,n_grid); 
n_it = 0;
bbeta = 1 / (1 + r);
change_price = 0 * eye(n_grid, n_grid); 
% price_1_chosen = 0 * eye(n_grid, n_grid); 
% price_2_chosen = 0 * eye(n_grid, n_grid); 


for jj = 1:n_grid

    diff = 1; 

    % The following is going to be inside a while 
    while diff > tol_level && n_it < max_n_it
        
        % Construct V_eroded by the inflation rate 
        V_eroded = 0 * eye(n_grid, n_grid); 
    
        for ii = 1:n_grid
    
           new_price_1 = z1_grid(ii,1) - g; 
           [min_value_1, index_min_1] = min(abs(new_price_1 - z1_grid)); 
    
           for kk = 1:n_grid 
    
                new_price_2 = z2_grid(kk,1) - g;
                [min_value_2, index_min_2] = min(abs(new_price_2 - z2_grid));
                
                V_eroded(ii,kk) = V_initial(index_min_1, index_min_2); 
            
           end
    
        end
        
        V_eroded(1,:) = -10^(5);
        V_eroded(:,1) = -10^(5);   
    
        % For loop 
    
        for ii = 1:n_grid
    
            %for jj = 1:n_grid 
    
                % Do nothing 
                new_price_1 = z1_grid(ii,1) - g;
                [~,new_price_1_index] = min(abs(z1_grid - new_price_1)); 
                new_price_2 = z2_grid(jj,1) - g;
                [~,new_price_2_index] = min(abs(z2_grid - new_price_2));
    
                V_doing_nothing = quadratic_profit(z1_grid(ii,1), z2_grid(jj,1), a, b, c) / 2 ...
                    + bbeta * V_initial(new_price_1_index, new_price_2_index);
                
                % Pay fixed cost and change prices 
                    % Grid search being careful
                    
                    V_changing_price = max(...
                        quadratic_profit(z1_grid, z2_grid(jj),a,b,c) / 2 + ...
                        bbeta * V_eroded(:,new_price_2_index)) ...
                        - K;
    
                    if V_changing_price > V_doing_nothing
    
                        change_price(ii,jj) = 1;
                        %price_1_chosen(ii,jj) = z1_grid(p1_chosen_index,1);
                        %price_2_chosen(ii,jj) = z2_grid(p2_chosen_index,1); 
    
                    else
    
                        change_price(ii,jj) = 0;
                        %price_1_chosen(ii,jj) = z1_grid(ii,1);
                        %price_2_chosen(ii,jj) = z2_grid(jj,1); 
    
                    end
    
                    V_new(ii,jj) = max(V_doing_nothing, V_changing_price); 
                  
            %end
    
        end
    
        
    
        diff = max(abs(V_new(:,jj) - V_initial(:,jj)));
        n_it = n_it + 1; 
        V_initial = V_new; 
    
    end

end

% surf(z1_grid, z2_grid, (1-bbeta)*V_new)
% plot(z1_grid, (1-bbeta)*V_new(:,25))

% Plots
figure(2)
plot(z1_grid, (1-bbeta) * V_new(:,floor(n_grid/4)), 'LineWidth',2)
hold on
plot(z1_grid, (1-bbeta) * V_new(:,floor(2*n_grid/4)), 'LineWidth',2)
plot(z1_grid, (1-bbeta) * V_new(:,floor(3*n_grid/4)), 'LineWidth',2)
xlabel('z1')
ylabel('Value Function')
title('Value Function Competitive Case with q = 0')
legend({strcat('z2= ', num2str(z1_grid(floor(n_grid/4)))), ...
    strcat('z2= ', num2str(z1_grid(floor(n_grid/2)))), ...
    strcat('z2= ', num2str(z1_grid(floor(3*n_grid/4))))}, ...
    'Location', 'northeast')
saveas(gcf, 'collusion_and_inflation_4_value_fun_competitive_0.png')
close(figure(2))


display((1-bbeta)*max(max(V_new)))




%%

% VFI when q=1 
% Step by step 
    % Suppose the state of the other firm is z2=1.5, let's do VFI of firm 1
    % only for that case -- I cannnot do it this way, because z2 is moving.
    % I have to do it together similar to the collusive case 

clearvars 
% Parameters
    % (real) Menu cost 
    K = 0.05;

    % Inflation rate
    g = 0.02; 

    % Interest rate 
    r = 0.03; 

    % Parameters of profit function
    a = 1; 
    c = 1/4; 
    b = 1/2; 

% VFI Competitive case with q = 1 -----------------------------------------

% Generate a grid
n_grid = 101;
tol_level = 10^(-2);
max_n_it = 1000; 
z1_grid = transpose(linspace(1,2,n_grid)); 
z2_grid = transpose(linspace(1,2,n_grid)); 
V_initial = 15 * ones(n_grid, n_grid); 
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
       [min_value_1, index_min_1] = min(abs(new_price_1 - z1_grid)); 

       for jj = 1:n_grid 

            new_price_2 = reaction_function(z1_grid(ii,1), z2_grid(jj,1), ...
                a, b, c, r, K) - g;
            [min_value_2, index_min_2] = min(abs(new_price_2 - z2_grid));
            
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

 

% Plots
figure(3)
plot(z1_grid, (1-bbeta) * V_new(:,floor(n_grid/4)), 'LineWidth',2)
hold on
plot(z1_grid, (1-bbeta) * V_new(:,floor(2*n_grid/4)), 'LineWidth',2)
plot(z1_grid, (1-bbeta) * V_new(:,floor(3*n_grid/4)), 'LineWidth',2)
xlabel('z1')
ylabel('Value Function')
title('Value Function Competitive Case with q = 1')
legend({strcat('z2= ', num2str(z1_grid(floor(n_grid/4)))), ...
    strcat('z2= ', num2str(z1_grid(floor(n_grid/2)))), ...
    strcat('z2= ', num2str(z1_grid(floor(3*n_grid/4))))}, ...
    'Location', 'northeast')
saveas(gcf, 'collusion_and_inflation_4_value_fun_competitive_1.png')
close(figure(3))


display((1-bbeta)*max(max(V_new)))






