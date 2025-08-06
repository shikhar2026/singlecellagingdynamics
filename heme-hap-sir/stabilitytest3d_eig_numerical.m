function nonlinear_system_stability()
    % Set random seed for reproducibility
    rng(1);
    
    % ===== SYSTEM PARAMETERS =====
    n = 10;
    k1 = 245.000; k4 = 0.22; KE = 3100.0;
    k7 = 2900.0; KM5 = 2900.0; KM1 = 280;
    KM3 = 2200; KM4 = 90; Stotal = 225.0;
    alpha = 0.35; beta = 0.95;
    k2 = 1.0; k3 = 0.004; k5 = 0.04;
    k6 = 0.004; k8 = 6.0; k9 = 1.0328;
    n1 = 2; n2p = 4; n3 = 2; n4 = 4; n5 = 6;

    % ===== SYSTEM DEFINITION =====
    odefun = @(t,X) [
        k1 * (X(3)^n2p / (X(3)^n2p + KE^n2p)) * ((1 - alpha) * X(1)^n1 + KM1) / (X(1)^n1 + KM1) + k2 - k3 * X(2);
        k4 * ((1 - beta) * X(2)^n3 + KM3) / (X(2)^n3 + KM3) * (X(1)^n4 / (X(1)^n4 + KM4^n4)) * (Stotal - X(1)) + k5 - k6 * X(1);
        k7 * (X(2)^n5 / (X(2)^n5 + KM5^n5)) + k8 - k9 * X(3)
    ];
    
    % ===== FIND FIXED POINTS =====
    fprintf('=== Searching for fixed points ===\n');
    
    % Search parameters
    Ntrials = 500;
    search_ranges = [0 600; 0 12000; 0 10000]; % [S; H; E] ranges
    fixed_points = [];
    
    options = optimoptions('fsolve','Display','off','TolFun',1e-10,'TolX',1e-10);
    
    for i = 1:Ntrials
        % Random initial guess
        X0 = search_ranges(:,1) + diff(search_ranges,1,2).*rand(3,1);
        
        % Find fixed point
        [X_fp, ~, exitflag] = fsolve(@(X) odefun(0,X), X0, options);
        
        % Validate solution
        if exitflag > 0 && all(isfinite(X_fp)) && all(X_fp >= 0)
            residual = norm(odefun(0,X_fp));
            if residual < 1e-6
                % Check uniqueness
                is_new = true;
                if ~isempty(fixed_points)
                    min_dist = min(vecnorm(fixed_points - X_fp', 2, 2));
                    is_new = (min_dist > 1e-3);
                end
                
                if is_new
                    fixed_points = [fixed_points; X_fp'];
                    fprintf('Found FP %d: [%.2f, %.2f, %.2f]\n', ...
                            size(fixed_points,1), X_fp');
                end
            end
        end
    end
    
    if isempty(fixed_points)
        fprintf('No fixed points found!\n');
        return;
    end
    
    % Remove duplicates
    fixed_points = uniquetol(fixed_points, 1e-3, 'ByRows', true);
    num_fps = size(fixed_points,1);
    
    % ===== STABILITY ANALYSIS =====
    fprintf('\n=== Stability Analysis ===\n');
    
    % Symbolic Jacobian
    syms S H E real
    F = [
        k1*(E^n2p/(E^n2p + KE^n2p))*((1-alpha)*S^n1 + KM1)/(S^n1 + KM1) + k2 - k3*H;
        k4*((1-beta)*H^n3 + KM3)/(H^n3 + KM3)*(S^n4/(S^n4 + KM4^n4))*(Stotal - S) + k5 - k6*S;
        k7*(H^n5/(H^n5 + KM5^n5)) + k8 - k9*E
    ];
    J_sym = jacobian(F, [S, H, E]);
    jacobian_fun = matlabFunction(J_sym, 'Vars', {S, H, E});
    
    % Numerical stability test parameters
    perturbation_size = 0.02; % 2% perturbation
    tspan = [0 2000]; % Long integration time
    colors = lines(5); % Different colors for each test trajectory
    
    for i = 1:num_fps
        fp = fixed_points(i,:)';
        fprintf('\nFixed Point %d/%d: [%.4f, %.4f, %.4f]\n', i, num_fps, fp);
        
        % Eigenvalue analysis
        J = jacobian_fun(fp(1), fp(2), fp(3));
        eigs = eig(J);
        fprintf('Eigenvalues:\n');
        disp(eigs');
        
        % Classify stability from eigenvalues
        re_eigs = real(eigs);
        if all(re_eigs < -1e-6)
            eig_stable = true;
            fprintf('Eigenvalue Prediction: Stable\n');
        else
            eig_stable = false;
            fprintf('Eigenvalue Prediction: Unstable\n');
        end
        
        % Numerical stability test
        fprintf('\nNumerical Stability Test:\n');
        num_stable = true;
        test_results = zeros(5,3);
        all_trajectories = cell(5,1);
        
        figure('Position', [100, 100, 1200, 500]);
        
        for test = 1:5
            % Create perturbed initial condition
            X0_pert = fp .* (1 + perturbation_size*(2*rand(3,1)-1));
            
            % Integrate system
            [t,X] = ode15s(@(t,X) odefun(t,X), tspan, X0_pert);
            final_state = X(end,:)';
            all_trajectories{test} = X;
            
            % Check convergence
            dist_to_fp = norm(final_state - fp);
            conv_ratio = dist_to_fp/norm(fp);
            
            test_results(test,:) = final_state';
            
            if conv_ratio > 0.05 % 5% tolerance
                num_stable = false;
            end
            
            fprintf('Test %d: Initial perturbation: [%.2f, %.2f, %.2f] -> Final state: [%.2f, %.2f, %.2f]\n', ...
                   test, X0_pert', final_state');
        end
        
        if num_stable
            stability_str = 'Stable';
        else
            stability_str = 'Unstable';
        end
        fprintf('Numerical Conclusion: %s\n', stability_str);
        
        % Combined stability assessment
        if eig_stable && num_stable
            fprintf('Combined Stability: Strongly Stable\n');
        elseif eig_stable && ~num_stable
            fprintf('Combined Stability: Theoretically Stable but Numerically Unstable\n');
        elseif ~eig_stable && num_stable
            fprintf('Combined Stability: Numerically Stable but Theoretically Unstable\n');
        else
            fprintf('Combined Stability: Confirmed Unstable\n');
        end
        
        % Visualization of test results with trajectories
        subplot(1,2,1);
        scatter3(fp(1), fp(2), fp(3), 100, 'r', 'filled'); hold on;
        
        % Plot each trajectory with different color
        for test = 1:5
            X_traj = all_trajectories{test};
            plot3(X_traj(:,1), X_traj(:,2), X_traj(:,3), 'Color', colors(test,:), 'LineWidth', 1.5);
            scatter3(X_traj(1,1), X_traj(1,2), X_traj(1,3), 50, colors(test,:), 'filled');
            scatter3(X_traj(end,1), X_traj(end,2), X_traj(end,3), 50, colors(test,:), 'filled', 'Marker', 'd');
        end
        
        xlabel('S'); ylabel('H'); zlabel('E');
        title(sprintf('FP %d Stability Test with Trajectories', i));
        legend_items = {'Fixed Point'};
        for test = 1:5
            legend_items = [legend_items, {sprintf('Test %d Trajectory', test), ...
                           sprintf('Test %d Start', test), sprintf('Test %d End', test)}];
        end
        legend(legend_items, 'Location', 'bestoutside');
        grid on; rotate3d on;
        
        subplot(1,2,2);
        plot(real(eigs), imag(eigs), 'bx', 'MarkerSize', 10); hold on;
        plot([0,0], ylim, 'k--');
        plot(xlim, [0,0], 'k--');
        xlabel('Real'); ylabel('Imaginary');
        title('Eigenvalues');
    end
end