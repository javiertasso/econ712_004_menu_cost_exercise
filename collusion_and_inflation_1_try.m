% -------------------------------------------------------------------------
% Collusion and Inflation 1 - First try 
% 2022
% Javiera Garcia & Javier Tasso
% S&W with two firms? Quadratic profit function and analysis following
% McMillan & Zinde-Walsh: Inflation and the timing of price changes
% -------------------------------------------------------------------------
clearvars
clc  
% cd 'C:\Users\Javier\Documents\Penn\Econ712-Macro-Financial-Markets\Collusion_and_Inflation'
cd 
addpath([cd '\collusion_and_inflation_1_try_functions'])
% -------------------------------------------------------------------------

% Set values of parameters 

    % Inflation rate 
    g = 0.02;

    % Interest rate 
    r= 0.03; 
    
    % (real) Menu cost
    bbeta = 0.05;
    
    % Number of firms 
    k = 2;

    % Monopoly price benchmark
    w_monop = log(1/2);
    profits_monop = instantaneous_profit(w_monop, w_monop + 0.01); 

% Try functions

    % Real price 
    real_price(g, 0.5, 2); 

    % Nominal price
    nominal_price(g, real_price(g, 0.5, 2),2); 
  
    % Instantaneous profit
    instantaneous_profit(-0.75, [-0.80; -0.72]); 
    instantaneous_profit(-0.75, [-0.72; -0.72]);
    instantaneous_profit(-0.75, [-0.75; -0.60]); 

    % Discounted profits between 0 and T when prices match 
    V_tilde(w_monop,w_monop,g,bbeta,r);
    V_tilde(w_monop-0.01,w_monop+0.01,g,bbeta,r);
    V_tilde(w_monop-0.05,w_monop+0.05,g,bbeta,r);

    % S&W threshold functions
    sheshinsky_weiss(exp(w_monop - 0.05), exp(w_monop + 0.05), ...
        r, g, bbeta);  

clear ans 

% What is the solution for Sheshinsky & Weiss? 
    % Assume they split demand in half
    % Try to find (s,S)
    [SW_thresholds,~] = fsolve(@(x) sheshinsky_weiss(x(1),x(2),r,g,bbeta), [exp(w_monop - 0.05),...
        exp(w_monop + 0.05)]);
    SW_thresholds_logs = log(SW_thresholds);  
 

% Now I want to find values a and b that make sense 
    % Problem, many trhesholds make sense here :( 
    % Guess some (a,b) such that a<b
        % 1. Check V_tilde(a,b)>0
        % 2. Check F_tilde(a) > 0
        % 3. Check V_tilde(a,b) > V_tilde(a,c) * disc FOR SOME c
            % I think this should be for any c too. 
        % 4. For any d<b-a the expression should be negative 
            % Here I can maximize the left hand side: if the maximum is
            % negative, then everything is negative. Use this to pin down
            % d=b-a? 

    % Choose (a,b)
        % 1. Fix V_tilde(a,b) a constant + F_tilde(a)>0
        % 2. Plot V_tilde(a,c) * e^(r(b-c)) with c in the horizontal axis I
        % want this function to be smaller than the constant for values of
        % c between a and b
        % 3. Plot the function in equation (18) with d in the horizontal
        % axis. I want this function to be smaller than the constant
        % whenever d<b-a

    % Suppose the following 
    a = SW_thresholds_logs(1,1);
    b = SW_thresholds_logs(1,2);
    %a = w_monop - 0.05; 
    %b = w_monop + 0.05;
    d_max = b-a; 
    c_min = a; 
    c_max = b; 
    const = V_tilde(a,b,g,bbeta,r); 
    F_tilde_a = instantaneous_profit(a,a); % I want this to be positive 

    % Start with condition involving c
    n_grid = 1001; 
    c_grid = transpose(linspace(c_min, c_max, n_grid)); 
    F_c = eye(n_grid,1); 
    
    for ii = 1:n_grid 

        F_c(ii,1) = exp(r*(b-c_grid(ii))) * V_tilde(a,c_grid(ii),g,bbeta,r); 

    end

    % Then condition involing d
    d_grid = transpose(linspace(0,d_max,n_grid)); 
    F_d = eye(n_grid, 1); 

    for ii = 1:n_grid

        F_d(ii,1) = V_monop(a-d_grid(ii),a,g,bbeta,r) - V_tilde(b-d_grid(ii), ...
            d_grid(ii),g,bbeta,r) + bbeta * (1 - exp(-(r/g) * d_grid(ii)));

    end





    
