% -------------------------------------------------------------------------
% Collusion and Inflation 5 - Summary 
% 2022
% Javiera Garcia & Javier Tasso
% Loop for different inflation rates 
% -------------------------------------------------------------------------

clearvars
clc 
% cd 'C:\Users\Javier\Documents\Penn\Econ712-Macro-Financial-Markets\Collusion_and_Inflation'
cd
addpath([cd '\collusion_and_inflation_5_summary_functions'])
try_functions = 0; 


%% Try Functions 

if try_functions == 1

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
    
    % VFI for the collusive case
    [max_V, z1_star, z2_star, z1_low, z1_high, z2_low, z2_high, ~, ...
        ~, ~, V_opt] = ...
        VFI_collusive_case(101, 10^(-2), 10000, a, b, c, r, K, g); 
    
    % Display Values
    display(max_V);
    display(z1_star);
    display(z2_star);
    display(z1_low);
    display(z1_high);
    display(z2_low);
    display(z2_high)
    display(V_opt); 
    
    % Clear 
    clear max_V z1_star z2_star z1_low z1_high z2_low z2_high V_opt
     
    % VFI for the competitive case to guess an initial value
    [max_V, z1_star, z2_star, z1_low, z1_high, z2_low, z2_high, V_fun, ...
        z1_grid, z2_grid, V_opt] = VFI_competitive_case_1(11, 10^(-3), ...
        10000, a, b, c, r, K, g, []);
    
    % Display Values
    display(max_V);
    display(z1_star);
    display(z2_star);
    display(z1_low);
    display(z1_high);
    display(z2_low);
    display(z2_high)
    display(V_opt); 
     
    % Clear 
    clear max_V z1_star z2_star z1_low z1_high z2_low z2_high V_opt
    
    % Generate a good initial guess 
    n_grid = 101;
    V_initial_guess = 0 * eye(n_grid, n_grid);
    z1_grid_b = transpose(linspace(1,2,n_grid));
    z2_grid_b = transpose(linspace(1,2,n_grid));
    
    % Fill up the initial guess using linear interpolation
    for jj = 1:n_grid
    
        for ii = 1:n_grid 
        
            % Get indices for z1
            [min_1, index_1] = min(abs(z1_grid_b(ii) - z1_grid)); 
            temp = abs(z1_grid_b(ii) - z1_grid); 
            temp(temp == min_1) = 10^5;
            [~, index_2] = min(temp); 
            index_1_min = min(index_1, index_2);
            index_1_max = index_1_min +1;
            clear min_1 index_1 temp index_2 
    
            % Get indices for z2
            [min_1, index_1] = min(abs(z2_grid_b(jj) - z2_grid)); 
            temp = abs(z2_grid_b(jj) - z2_grid); 
            temp(temp == min_1) = 10^5;
            [~, index_2] = min(temp); 
            index_2_min = min(index_1, index_2);
            index_2_max = index_2_min +1;
            clear min_1 index_1 temp index_2 
    
            V_initial_guess(ii,jj) = V_fun(index_1_min, index_2_min) + ...
                (V_fun(index_1_max, index_2_min) - V_fun(index_1_min, index_2_min)) / ...
                (z1_grid(index_1_max) - z1_grid(index_1_min)) * ...
                (z1_grid_b(ii) - z1_grid(index_1_min)) + ...
                (V_fun(index_1_min, index_2_max) - V_fun(index_1_min, index_2_min)) / ...
                (z2_grid(index_2_max) - z2_grid(index_2_min)) * ...
                (z2_grid_b(jj) - z2_grid(index_2_min));
        
        end
    
    end
     
    clear ii jj z1_grid_b z2_grid_b 
    
    % Try function with initial guess
    [max_V, z1_star, z2_star, z1_low, z1_high, z2_low, z2_high, V_fun, ...
        z1_grid, z2_grid, V_opt] = VFI_competitive_case_1(n_grid, 10^(-2), ...
        10000, a, b, c, r, K, g, V_initial_guess);
    
    % Display Values
    display(max_V);
    display(z1_star);
    display(z2_star);
    display(z1_low);
    display(z1_high);
    display(z2_low);
    display(z2_high)
    display(V_opt); 
    display(V_fun(z1_grid == z1_high, z2_grid == z2_high))
    display(V_fun(z1_grid == z1_low, z2_grid == z2_low))
    
    clear index_1_max index_1_min index_2_max index_2_min max_V n_grid V_fun ...
        V_initial_guess V_opt z1_grid z1_high z1_low z1_star z2_grid z2_high ...
        z2_low z2_star 
    
    % VFI_competitive_case_0
    [max_V, z1_star, z2_star, z1_low, z1_high, z2_low, z2_high, V_fun, ...
        z1_grid, z2_grid, V_opt] = VFI_competitive_case_0(101, 10^(-2), ...
        10000, a, b, c, r, K, g);
    
    % Display Values
    display(max_V);
    display(z1_star);
    display(z2_star);
    display(z1_low);
    display(z1_high);
    display(z2_low);
    display(z2_high)
    display(V_opt); 
    display(V_fun(z1_grid == z1_high, z2_grid == z2_high))
    display(V_fun(z1_grid == z1_low, z2_grid == z2_low))
    
    clearvars -except a b c g K r

end
%% Loop for some inflation rates 

% Given the above functions I want to do the following 
   % 1) For an inflation rate, get the value functions in each case
   % 2) Calculate the value of deviating 
   % 3) Figure out the delta that sustains collusion (careful, do I use
   % each period or do I aggregate periods? I think that the former is more
   % straightforward but it involves that the equation over delta will
   % depend on the time span of every inaction set.    

% Parameters - These ones are going to stay the same 
    % (real) Menu cost 
    K = 0.05;

    % Interest rate 
    r = 0.03; 
    bbeta = 1 / (1+r); 

    % Parameters of profit function
    a = 1; 
    c = 1/4; 
    b = 1/2; 

    % Parameters to initialize the algortihm 
    n_grid = 101;
    tol_level = 10^(-3); 
    max_n_it = 1000; 

    % Inflation Vector
    max_k = 10; 
    g_vector = transpose(linspace(0.01,0.05,max_k)); 
    

    % Store things 
    delta = zeros(max_k,1);
    V_colu = zeros(max_k,1);
    V_comp = zeros(max_k,1);
    z1_deviation = zeros(max_k,1);
    length_cycle = zeros(max_k,1);
    discounted_profits_devi = zeros(max_k,1); 
    discounted_profits_colu = zeros(max_k,1); 
    z1_low_colu = zeros(max_k,1); 
    z1_high_colu = zeros(max_k,1); 
    z1_low_comp_1 = zeros(max_k,1); 
    z1_high_comp_1 = zeros(max_k,1); 



for kk = 1:max_k
    
    % Display
    disp(kk); 

    % Inflation rate - Loop over different inflation rates 
    g = g_vector(kk); 
    
    % Get the value function when there is collusion 
    [~, ~, ~, z1_low_colu(kk), z1_high_colu(kk), ~, z2_high_colu, V_fun_colu, ...
        z1_grid_colu, z2_grid_colu, ~] = ...
        VFI_collusive_case(n_grid, tol_level, max_n_it, a, b, c, r, K, g); 
    V_fun_colu = V_fun_colu / 2; 
    
    % Get the value function when there is competition and the other firm
    % chooses q=0
    [~, ~, ~, ~, ~, ~, ~, ~, ...
        ~, ~, ~] = ...
        VFI_competitive_case_0(n_grid, tol_level, max_n_it, a, b, c, r, K, g); 
    
    % Get the value function when there is competition and the other firm
    % chooses q=1
    
        % Step 1: get an initial guess ----------------------------------------
        [~, ~, ~, ~, ~, ~, ~, V_fun, ...
        z1_grid, z2_grid, ~] = VFI_competitive_case_1(11, 10^(-3), ...
        max_n_it, a, b, c, r, K, g, []);
    
        % Generate a good initial guess 
        V_initial_guess = 0 * eye(n_grid, n_grid);
        z1_grid_b = transpose(linspace(1,2,n_grid));
        z2_grid_b = transpose(linspace(1,2,n_grid));
        
        % Fill up the initial guess using linear interpolation
        for jj = 1:n_grid
        
            for ii = 1:n_grid 
            
                % Get indices for z1
                [min_1, index_1] = min(abs(z1_grid_b(ii) - z1_grid)); 
                temp = abs(z1_grid_b(ii) - z1_grid); 
                temp(temp == min_1) = 10^5;
                [~, index_2] = min(temp); 
                index_1_min = min(index_1, index_2);
                index_1_max = index_1_min +1;
                clear min_1 index_1 temp index_2 
        
                % Get indices for z2
                [min_1, index_1] = min(abs(z2_grid_b(jj) - z2_grid)); 
                temp = abs(z2_grid_b(jj) - z2_grid); 
                temp(temp == min_1) = 10^5;
                [~, index_2] = min(temp); 
                index_2_min = min(index_1, index_2);
                index_2_max = index_2_min +1;
                clear min_1 index_1 temp index_2 
        
                V_initial_guess(ii,jj) = V_fun(index_1_min, index_2_min) + ...
                    (V_fun(index_1_max, index_2_min) - V_fun(index_1_min, index_2_min)) / ...
                    (z1_grid(index_1_max) - z1_grid(index_1_min)) * ...
                    (z1_grid_b(ii) - z1_grid(index_1_min)) + ...
                    (V_fun(index_1_min, index_2_max) - V_fun(index_1_min, index_2_min)) / ...
                    (z2_grid(index_2_max) - z2_grid(index_2_min)) * ...
                    (z2_grid_b(jj) - z2_grid(index_2_min));
            
            end
        
        end
    
        clear ii jj index_1_max index_1_min index_2_max index_2_min V_fun z1_grid ...
            z1_grid_b z2_grid z2_grid_b 
    
        % Step 2: do VFI with the initial guess -------------------------------
        [~, ~, ~, z1_low_comp_1(kk), z1_high_comp_1(kk), ~, z2_high_comp_1, V_fun_comp_1, ...
        ~, ~, ~] = VFI_competitive_case_1(n_grid, tol_level, ...
        max_n_it, a, b, c, r, K, g, V_initial_guess);
    
        clear V_initial_guess 
    
    % Check wether collusion and competition set different prices
    % if z1_high_colu(kk) == z1_high_comp_1
    
        % disp('Collusion and competition set the same price')
    
    % else
    
        % disp('Collusion and competition set different prices') 
    
    % end
    
    % Number of periods the price is going to stay fixed 
    length_cycle(kk) = round((z1_high_comp_1(kk) - z1_low_comp_1(kk)) / g); 
    z1_grid = z1_grid_colu; 
    z2_grid = z2_grid_colu; 
    
    % Amount of profits in that cycle 
    profit_per_period_colu = 0 * eye(length_cycle(kk), 1); 
    profit_per_period_devi = 0 * eye(length_cycle(kk), 1); 
    discount_vector = 0 * eye(1, length_cycle(kk)); 
    
    for tt = 1:length_cycle(kk) 
    
        profit_per_period_colu(tt) = quadratic_profit(z1_high_colu(kk) - (tt-1) * g, ...
            z2_high_colu - (tt-1) * g, a, b, c) / 2; 
        discount_vector(tt) = bbeta^(tt-1); 
    
    end
    
    % Profits for staying in that collusive agreement 
    discounted_profits_colu(kk) = discount_vector * profit_per_period_colu; 
    
    clear tt z1_grid_colu z2_grid_colu 
        
    % Profits for deviating during that period 
    profit_per_period_devi_temp = 0 * eye(n_grid, 1); 
    
    for ii = 1:n_grid 
    
        for tt = 1:length_cycle(kk)
    
            profit_per_period_devi(tt) = quadratic_profit(z1_grid(ii) - (tt-1) * g, ...
                z2_high_colu - (tt-1) * g, a, b, c) / 2;
    
        end
    
         profit_per_period_devi_temp(ii) = discount_vector * profit_per_period_devi; 
    
    end
    
    [discounted_profits_devi(kk), index_devi] = max(profit_per_period_devi_temp); 
    
    z1_deviation(kk) = z1_grid(index_devi); 
    
    clear ii tt profit_per_period_devi_temp index_devi
    
    % Outputs
    delta(kk) = ((discounted_profits_devi(kk) - discounted_profits_colu(kk)) / ... 
        (V_fun_colu(z1_grid == z1_high_colu(kk), z2_grid == z2_high_colu) / (1-bbeta) - ...
        V_fun_comp_1(z1_grid == z1_high_comp_1(kk), z2_grid == z2_high_comp_1) / (1-bbeta)))^...
        (1/length_cycle(kk));
    
    V_colu(kk) = V_fun_colu(z1_grid == z1_high_colu(kk), z2_grid == z2_high_colu) / (1-bbeta); 
    V_comp(kk) = V_fun_comp_1(z1_grid == z1_high_comp_1(kk), z2_grid == z2_high_comp_1) / (1-bbeta); 

end

%% Store results 

results_inflation_loop = [g_vector, length_cycle, discounted_profits_devi, ...
    discounted_profits_colu, delta, V_colu, V_comp, z1_deviation, z1_low_colu, ...
    z1_high_colu, z1_low_comp_1, z1_high_comp_1]; 
colNames = {'Inflation', 'Length', 'Profit_Dev', 'Profit_Col', 'Delta', ...
    'V_colu', 'V_comp', 'z_dev', 'z_low_colu', 'z_high_colu', 'z_low_comp', ...
    'z_high_comp'};
sTable = array2table(results_inflation_loop, 'VariableNames', colNames); 
writetable(sTable, 'collusion_and_inflation_5_results.csv')

%% Plot results 
clearvars
r = 0.03;
bbeta = 1 / (1+r);
clear r;

% Import results 

    % Set up the Import Options and import the data
    opts = delimitedTextImportOptions("NumVariables", 12);
    
    % Specify range and delimiter
    opts.DataLines = [2, Inf];
    opts.Delimiter = ",";
    
    % Specify column names and types
    opts.VariableNames = ["Inflation", "Length", "Profit_Dev", "Profit_Col", "Delta", "V_colu", "V_comp", "z_dev", "z_low_colu", "z_high_colu", "z_low_comp", "z_high_comp"];
    opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];
    
    % Specify file level properties
    opts.ExtraColumnsRule = "ignore";
    opts.EmptyLineRule = "read";
    
    % Import the data
    results = readtable("C:\Users\Javier\Documents\Penn\Econ712-Macro-Financial-Markets\Collusion_and_Inflation\collusion_and_inflation_5_results.csv", opts);
    
    % Clear temporary variables
    clear opts

% Change directory to save plots 
% cd 'C:\Users\Javier\Documents\Penn\Econ712-Macro-Financial-Markets\Collusion_and_Inflation'
% cd 'C:\Users\Javier\Dropbox\Aplicaciones\Overleaf\Collusion_Inflation'


% Plots
figure(1)
plot(results.Inflation, results.z_high_colu, 'LineWidth',2,'Color', 'b')
hold on
plot(results.Inflation, results.z_high_comp, 'LineWidth',2,'Color', 'r')
plot(results.Inflation, results.z_dev, 'LineWidth',2,'Color', 'g')
plot(results.Inflation, results.z_low_colu, 'LineWidth',2,'Color', 'b')
plot(results.Inflation, results.z_low_comp, 'LineWidth',2,'Color', 'r')
legend('Collusion','Competition', 'Deviation', 'Location','best')
title('Inaction Sets for Different Inflation Rates')
xlabel('Inflation')
ylabel('Thresholds')
saveas(gcf, 'collusion_and_inflation_5_thresholds_colu_comp.png')
close(figure(1))

figure(2)
plot(results.Inflation, results.Length, 'LineWidth',2,'Color', 'b')
hold on
title('Time Between Price Changes')
xlabel('Inflation')
ylabel('Periods')
saveas(gcf, 'collusion_and_inflation_5_time_inaction.png')
close(figure(2))

figure(3)
plot(results.Inflation, (1-bbeta) * results.V_colu, 'LineWidth',2,'Color', 'b')
hold on
plot(results.Inflation, (1-bbeta) * results.V_comp, 'LineWidth',2,'Color', 'r')
legend('Collusion','Competition', 'Location','best')
title('Value of Each Alternative')
xlabel('Inflation')
ylabel('Value')
saveas(gcf, 'collusion_and_inflation_5_value.png')
close(figure(3))
 
figure(4)
plot(results.Inflation, results.Delta, 'LineWidth',2,'Color', 'b')
hold on
title('Discount Factor for Different Inflation Rates')
xlabel('Inflation')
ylabel('Discount Factor')
saveas(gcf, 'collusion_and_inflation_5_delta.png')
close(figure(4))

figure(5)
plot(results.Inflation, results.Profit_Col ./ results.Length, 'LineWidth',2,'Color', 'b')
hold on
plot(results.Inflation, results.Profit_Dev ./ results.Length, 'LineWidth',2,'Color', 'g')
legend('Collusion','Deviation', 'Location','best')
title('Average Profit Over the Period')
xlabel('Inflation')
ylabel('Average Profits')
saveas(gcf, 'collusion_and_inflation_5_mean_profit.png')
close(figure(5))

clearvars
exit 





