% Provided data
E_data = [2,6,17,38,75,150,300];   % heme (independent for H vs E)
H_data = [0.5,0.5,3.2,16,39,67,124]; % hap  (dependent for H vs E)

% --- Fit H as a function of E: H = Vmax1 * E^n2 / (Km2^n2 + E^n2) ---
hillFun_HE = @(params, E) params(1) * (E.^params(3)) ./ (params(2).^params(3) + E.^params(3));
initParams_HE = [max(H_data), median(E_data), 2]; % [Vmax1, Km2, n2]
lb_HE = [0, 0, 0.5];
ub_HE = [2*max(H_data), 2*max(E_data), 10];

params_HE = lsqcurvefit(hillFun_HE, initParams_HE, E_data, H_data, lb_HE, ub_HE);

Vmax1 = params_HE(1);
Km2   = params_HE(2);
n2    = params_HE(3);

fprintf('Fit for H = Vmax1 * E^{n2} / (Km2^{n2} + E^{n2}):\n');
fprintf('  Vmax1 = %.4f\n  Km2 = %.4f\n  n2 = %.4f\n', Vmax1, Km2, n2);

% --- Fit E as a function of H: E = Vmax2 * H^n5 / (Km5^n5 + H^n5) ---
hillFun_EH = @(params, H) params(1) * (H.^params(3)) ./ (params(2).^params(3) + H.^params(3));
initParams_EH = [max(E_data), median(H_data), 2]; % [Vmax2, Km5, n5]
lb_EH = [0, 0, 0.5];
ub_EH = [2*max(E_data), 2*max(H_data), 10];

params_EH = lsqcurvefit(hillFun_EH, initParams_EH, H_data, E_data, lb_EH, ub_EH);

Vmax2 = params_EH(1);
Km5   = params_EH(2);
n5    = params_EH(3);

fprintf('Fit for E = Vmax2 * H^{n5} / (Km5^{n5} + H^{n5}):\n');
fprintf('  Vmax2 = %.4f\n  Km5 = %.4f\n  n5 = %.4f\n', Vmax2, Km5, n5);

% --- Plotting both fits ---
E_fit = linspace(min(E_data), max(E_data), 200);
H_fit = linspace(min(H_data), max(H_data), 200);

figure;
subplot(1,2,1);
plot(E_data, H_data, 'ko', 'MarkerFaceColor','k'); hold on;
plot(E_fit, hillFun_HE(params_HE, E_fit), 'r-', 'LineWidth',2);
xlabel('E (heme)'); ylabel('H (hap)');
title('H = Hill(E)');
legend('Data','Hill Fit');

subplot(1,2,2);
plot(H_data, E_data, 'ko', 'MarkerFaceColor','k'); hold on;
plot(H_fit, hillFun_EH(params_EH, H_fit), 'b-', 'LineWidth',2);
xlabel('H (hap)'); ylabel('E (heme)');
title('E = Hill(H)');
legend('Data','Hill Fit');
% After fitting, for H as a function of E:
Vmax1 = params_HE(1);   % Vmax for H = Hill(E)

% For E as a function of H:
Vmax2 = params_EH(1);   % Vmax for E = Hill(H)

fprintf('Vmax for H = Hill(E): %.4f\n', Vmax1);
fprintf('Vmax for E = Hill(H): %.4f\n', Vmax2);



% Known parameters from previous fit
k7 = 600;      % Vmax2 from your previous Hill fit
n5 = 1.6158;      % Hill coefficient from previous fit
Km5 = 126.6173;     % Km from previous fit

% Model function for E as a function of H, k8, and k9
modelFun = @(params, H) (k7 * (H.^n5 ./ (H.^n5 + Km5.^n5)) + params(1)) ./ params(2);
% params(1) = k8, params(2) = k9

% Initial guesses for [k8, k9]
initParams = [1, 1];

% Fit k8 and k9 using lsqcurvefit
paramsFit = lsqcurvefit( ...
    modelFun, ...
    initParams, ...
    H_data, ...
    E_data ...
);

k8 = paramsFit(1);
k9 = paramsFit(2);

fprintf('Fitted steady-state parameters:\n');
fprintf('k8 = %.4f\n', k8);
fprintf('k9 = %.4f\n', k9);

% Optional: Plot fit
H_fit = linspace(min(H_data), max(H_data), 200);
E_fit = modelFun(paramsFit, H_fit);

figure;
plot(H_data, E_data, 'ko', 'MarkerFaceColor', 'k'); hold on;
plot(H_fit, E_fit, 'r-', 'LineWidth', 2);
xlabel('H (HAP)'); ylabel('E (heme)');
title('Steady-State Fit: E as function of H');
legend('Data', 'Model Fit');
grid on;
