rng(1);
delete(gcp('nocreate'));

parpool(10)
% --- Target Fixed Points (H, S) ---
target_HS = [
    250.107, 187.207;
    250.119, 37.7391;
    250.120, 10.1732;
    3568.41, 10.0521;
    3610.95, 66.9847;
    3713.01, 126.417;
    9919.71, 10.0156
];

% --- Parameter Ranges (broadened) ---
k1_range  = 245:5:255;
k4_range  = 0.18:0.02:0.22;
KE_range  = 2900:100:3100;
k7_range  = 2900:100:3100;
KM5_range = 2900:100:3100;
n_range   = 6:2:10; % Hill coefficients independently

% --- Other Parameters (fixed) ---
alpha = 0.35; KM1 = 280; k2 = 1.0; k3 = 0.004;
beta = 0.95; KM3 = 2200; KM4 = 90;
Stotal = 225.0; k5 = 0.04; k6 = 0.004;
k8 = 6.0; k9 = 1.0328;

tol_resid = 1e-2;
tol_unique = 1e-1;

Ntrials = 100; % Increased trials for better fixed point finding

% --- For recording best parameter set ---
best_total_dist = Inf;
best_results = struct('params', [], 'nvals', [], 'total_dist', [], 'pairs', [], 'fixed_points', []);

% Start scanning parameters
for n1 = n_range
for n2p = n_range
for n3 = n_range
for n4 = n_range
for n5 = n_range
    for k1 = k1_range
        for k4 = k4_range
            for KE = KE_range
                for k7 = k7_range
                    for KM5 = KM5_range
                        % Define ODE function for current parameters
                        odefun = @(X) [
                            k1 * (X(3)^n2p / (X(3)^n2p + KE^n2p)) * ((1 - alpha) * X(1)^n1 + KM1) / (X(1)^n1 + KM1) + k2 - k3 * X(2);
                            k4 * ((1 - beta) * X(2)^n3 + KM3) / (X(2)^n3 + KM3) * (X(1)^n4 / (X(1)^n4 + KM4^n4)) * (Stotal - X(1)) + k5 - k6 * X(1);
                            k7 * (X(2)^n5 / (X(2)^n5 + KM5^n5)) + k8 - k9 * X(3)
                        ];

                        fixed_points = [];
                        Xmin = [0, 0, 0]; Xmax = [600, 12000, 10000];
                        options = optimoptions('fsolve','Display','off','TolFun',1e-10,'TolX',1e-10);

                        % Find fixed points with multiple random initial guesses
                        for i = 1:Ntrials
                            X0 = Xmin + (Xmax - Xmin).*rand(1,3);
                            [X_fp, fval, exitflag] = fsolve(odefun, X0, options);
                            if exitflag > 0 && all(isreal(X_fp)) && all(X_fp >= 0) && all(abs(fval) < tol_resid)
                                X_fp = real(X_fp);
                                if isempty(fixed_points)
                                    fixed_points = X_fp;
                                else
                                    diffs = abs(fixed_points - X_fp);
                                    dists = sum(diffs, 2);
                                    if all(dists > tol_unique)
                                        fixed_points = [fixed_points; X_fp];
                                    end
                                end
                            end
                        end

                        if isempty(fixed_points)
                            % No fixed points found, skip
                            continue;
                        end

                        % Remove duplicates and round
                        fixed_points = round(fixed_points, 6, 'significant');
                        [~, ia, ~] = unique(fixed_points, 'rows');
                        fixed_points = fixed_points(ia, :);

                        % Check if enough fixed points for assignment
                        if size(fixed_points,1) < size(target_HS,1)
                            % Not enough fixed points, skip
                            continue;
                        end

                        % Compute distance matrix (targets x fixed points)
                        found_HS = fixed_points(:, [2,1]); % [H, S]
                        D = pdist2(target_HS, found_HS);

                        % Solve assignment problem using Hungarian algorithm
                        if exist('matchpairs','file')
                            [pairs, total_dist] = matchpairs(D, 1e6);
                        else
                            % Greedy fallback assignment
                            pairs = zeros(size(target_HS,1),2);
                            used = false(1,size(found_HS,1));
                            total_dist = 0;
                            for ti = 1:size(target_HS,1)
                                [minval, idx] = min(D(ti,:) + used*1e6);
                                pairs(ti,:) = [ti, idx];
                                total_dist = total_dist + minval;
                                used(idx) = true;
                            end
                        end

                        % Print progress info
                        fprintf('\nn1=%d n2p=%d n3=%d n4=%d n5=%d, k1=%.1f, k4=%.2f, KE=%.1f, k7=%.1f, KM5=%.1f | Fixed points: %d | Sum dist: %.4f\n', ...
                            n1, n2p, n3, n4, n5, k1, k4, KE, k7, KM5, size(fixed_points,1), total_dist);
                        fprintf('Target (H,S) -> Fixed Point (S,H,E) | Distance\n');
                        for k = 1:size(pairs,1)
                            ti = pairs(k,1);
                            fi = pairs(k,2);
                            fp = fixed_points(fi,:);
                            dist = norm(target_HS(ti,:) - found_HS(fi,:));
                            fprintf('(%.3f, %.3f) -> (%.3f, %.3f, %.3f) | %.4f\n', ...
                                target_HS(ti,1), target_HS(ti,2), fp(1), fp(2), fp(3), dist);
                        end

                        % Update best result if improved
                        if total_dist < best_total_dist
                            best_total_dist = total_dist;
                            best_results.params = [k1, k4, KE, k7, KM5];
                            best_results.nvals = [n1, n2p, n3, n4, n5];
                            best_results.total_dist = total_dist;
                            best_results.pairs = pairs;
                            best_results.fixed_points = fixed_points;
                            best_results.found_HS = found_HS;
                        end
                    end
                end
            end
        end
    end
end
end
end
end
end

% --- Final report ---
fprintf('\n\n=== OPTIMIZED PARAMETER SET ===\n');
if isinf(best_total_dist)
    fprintf('No suitable parameter set found. Try increasing Ntrials or broadening parameter ranges.\n');
else
    p = best_results.params;
    nvals = best_results.nvals;
    fprintf('n1=%d, n2p=%d, n3=%d, n4=%d, n5=%d\n', nvals(1), nvals(2), nvals(3), nvals(4), nvals(5));
    fprintf('k1=%.1f, k4=%.2f, KE=%.1f, k7=%.1f, KM5=%.1f\n', p(1), p(2), p(3), p(4), p(5));
    fprintf('Minimum total Euclidean distance: %.4f\n', best_total_dist);
    fprintf('Target (H,S) -> Assigned Fixed Point (S,H,E) | Distance\n');
    for k = 1:size(best_results.pairs,1)
        ti = best_results.pairs(k,1);
        fi = best_results.pairs(k,2);
        fp = best_results.fixed_points(fi,:);
        dist = norm(target_HS(ti,:) - best_results.found_HS(fi,:));
        fprintf('(%.3f, %.3f) -> (%.3f, %.3f, %.3f) | %.4f\n', ...
            target_HS(ti,1), target_HS(ti,2), fp(1), fp(2), fp(3), dist);
    end
    fprintf('\nAll fixed points for this parameter set ([S, H, E]):\n');
    disp(best_results.fixed_points);
end
