function fixed_points = find_all_fixed_points()
    % Initialize parameters
    params = struct();
    params.k1 = 245.000;  params.k4 = 0.22;    params.KE = 3100.0;
    params.k7 = 2900.0;   params.KM5 = 2900.0;  params.KM1 = 280;
    params.KM3 = 2200;    params.KM4 = 90;      params.Stotal = 225.0;
    params.alpha = 0.35;  params.beta = 0.95;
    params.k2 = 1.0;      params.k3 = 0.004;    params.k5 = 0.04;
    params.k6 = 0.004;    params.k8 = 6.0;      params.k9 = 1.0328;
    params.n = 10;        params.n1 = params.n; params.n2p = params.n; 
    params.n3 = 6;        params.n4 = 6;        params.n5 = 8;

    % Define the system equations
    syms S H E real
    X = [S; H; E];
    
    f1 = params.k1*(E^params.n2p/(E^params.n2p + params.KE^params.n2p))*...
         ((1-params.alpha)*S^params.n1 + params.KM1)/(S^params.n1 + params.KM1) + ...
         params.k2 - params.k3*H;
    
    f2 = params.k4*((1-params.beta)*H^params.n3 + params.KM3)/...
         (H^params.n3 + params.KM3)*(S^params.n4/(S^params.n4 + params.KM4^params.n4))*...
         (params.Stotal - S) + params.k5 - params.k6*S;
    
    f3 = params.k7*(H^params.n5/(H^params.n5 + params.KM5^params.n5)) + ...
         params.k8 - params.k9*E;
    
    F = [f1; f2; f3];
    
    % Convert to numerical function
    odefun = matlabFunction(F, 'Vars', {X});
    
    % Search settings
    search_ranges = [0 600; 0 12000; 0 10000];  % [S_min S_max; H_min H_max; E_min E_max]
    N_trials = 500;          % Number of random starting points
    tol = 1e-4;              % Tolerance for considering a solution valid
    min_distance = 1;        % Minimum distance between distinct fixed points
    
    % Multi-start search
    fixed_points = [];
    options = optimoptions('fsolve', 'Display', 'off', 'FunctionTolerance', 1e-8);
    
    for i = 1:N_trials
        % Generate random initial guess within bounds
        X0 = search_ranges(:,1) + rand(3,1).*(search_ranges(:,2)-search_ranges(:,1));
        
        try
            % Find fixed point
            [X_fp, ~, exitflag] = fsolve(odefun, X0, options);
            
            % Check if solution is valid
            if exitflag > 0 && all(X_fp >= -1e-3) && norm(odefun(X_fp)) < tol
                % Round to 4 decimal places for comparison
                X_rounded = round(X_fp', 4);
                
                % Check if this point is new
                if isempty(fixed_points) || ...
                   min(vecnorm(fixed_points - X_rounded, 2, 2)) > min_distance
                    fixed_points = [fixed_points; X_rounded];
                end
            end
        catch
            continue
        end
    end
    
    % Display results
    fprintf('Found %d fixed points:\n', size(fixed_points,1));
    disp(fixed_points);
    
    % Visualize in 3D space
    figure;
    scatter3(fixed_points(:,1), fixed_points(:,2), fixed_points(:,3), 100, 'filled');
    xlabel('S'); ylabel('H'); zlabel('E');
    title(sprintf('Found %d Fixed Points', size(fixed_points,1)));
    grid on; rotate3d on;
end