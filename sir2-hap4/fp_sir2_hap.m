% --- Place your parameter values here ---
k1    = 40.0;      % e.g., 40
alpha = 0.35;      % e.g., 0.35
KM1   = 280.0;      % e.g., 280
n1    = 2.0;      % e.g., 2
k2    = 1.0;      % e.g., 1.0
k3    = 0.004;      % e.g., 0.004
n2    = 4.0;      % e.g., 4
KM2   = 4250.0;      % e.g., 4250
k4    = 0.02;      % e.g., 0.02
beta  = 0.95;      % e.g., 0.95
KM3   = 2200.0;      % e.g., 2200
n3    = 2.0;      % e.g., 2
KM4   = 90.0;      % e.g., 90
n4    = 4.0;      % e.g., 4
Stotal= 225.0;      % e.g., 40
k5    = 0.04;      % e.g., 0.04
k6    = 0.004;      % e.g., 0.004

% --- System of equations as a function handle ---
odefun = @(X) [
    k1 * ( (1 - alpha) * X(2)^n1 + KM1^n1 ) / ( X(2)^n1 + KM1^n1 ) * ...
        X(1)^n2 / ( X(1)^n2 + KM2^n2 ) + k2 - k3*X(1);
    k4 * ( (1 - beta) * X(1)^n3 + KM3^n3 ) / ( X(1)^n3 + KM3^n3 ) * ...
        X(2)^n4 / ( X(2)^n4 + KM4^n4 ) * ( Stotal - X(2) ) + k5 - k6*X(2)
];

% --- Settings for root search ---
Ntrials = 2000; % Number of random initial guesses (increase if needed)
Xmin = [0, 0];  % Lower bounds for [H, S] (positive only)
Xmax = [20000,1000]; % Upper bounds for [H, S] (adjust as needed)
tol_resid = 1e-6;      % Tolerance for residuals
tol_unique = 1e-4;     % Tolerance for uniqueness

fixed_points = [];
options = optimoptions('fsolve','Display','off','TolFun',1e-12,'TolX',1e-12);

for i = 1:Ntrials
    % Generate random initial guess in the box
    X0 = Xmin + (Xmax - Xmin).*rand(1,2);
    [X_fp, fval, exitflag] = fsolve(odefun, X0, options);
    % Accept if converged, real, positive, and small residuals
    if exitflag > 0 && all(isreal(X_fp)) && all(X_fp >= 0) && all(abs(fval) < tol_resid)
        X_fp = real(X_fp); % Remove any tiny imaginary part
        % Uniqueness check (rounded to 6 decimals)
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

% --- Remove any remaining near-duplicates (final pass) ---
fixed_points = round(fixed_points, 6, 'significant'); % Round for uniqueness
[~, ia, ~] = unique(fixed_points, 'rows');
fixed_points = fixed_points(ia, :);

% --- Print the unique real, positive fixed points ---
fprintf('\nUnique real, positive fixed points found:\n');
for i = 1:size(fixed_points,1)
    fprintf('Fixed point %d: H = %.6f, S = %.6f\n', i, fixed_points(i,1), fixed_points(i,2));
end

fprintf('\nTotal unique real, positive fixed points found: %d\n', size(fixed_points,1));
