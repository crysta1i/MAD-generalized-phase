function [params, omega, K_eff, fitinfo] = fit_kuramoto_structured(allses_angle_cts, ch_names, trial, fs) % chan_pos
% FIT_KURAMOTO_STRUCTURED
% Fit a Kuramoto model with structured, time-varying coupling:
%
%   dtheta_i/dt = omega_i + A(t) * sum_j exp(-d_ij/l) * sin(theta_j - theta_i)
%
% where
%   A(t) = A0 * exp(-(t - T0)^2 / (2*sigma^2))
%
% Inputs:
%   allses_angle_cts : (n_timepoints x n_contacts x n_trials) phase data
%   ch_names         : cell array of channel names, length N
%   trial            : trial index
%   fs               : sampling frequency in Hz
%   chan_pos         : N x D matrix of channel coordinates (D=2 or 3)
%
% Outputs:
%   params           : struct with fitted global parameters
%   omega            : N x 1 fitted natural frequencies (rad/s)
%   K_eff            : N x N effective coupling matrix at peak A(t)
%                      K_eff(i,j) = A0 * exp(-d_ij/l)
%   fitinfo          : optimization output struct
%

    N = numel(ch_names);

    % Frequency bounds for omega in rad/s
    f_lo = 5;
    f_hi = 40;
    omega_lo = 2*pi*f_lo;
    omega_hi = 2*pi*f_hi;

    %% 1) Load phase angles
    fprintf('Loading phase angles for %d channels, trial %d ...\n', N, trial);

    theta_cell = cell(N, 1);
    for i = 1:N
        theta_cell{i} = allses_angle_cts(:, i, trial);
    end

    for i = 1:N
        theta_cell{i} = theta_cell{i}(:)';  % row vector
    end

    T_lengths = cellfun(@numel, theta_cell);
    if ~all(T_lengths == T_lengths(1))
        error('Channels have different numbers of time points.');
    end

    T = T_lengths(1);
    theta = vertcat(theta_cell{:});   % N x T

    %% 2) Unwrap and estimate derivatives
    theta_u = unwrap(theta, [], 2);
    dTheta = diff(theta_u, 1, 2) * fs;     % N x (T-1)

    % Midpoint phases for RHS
    theta_mid = (theta_u(:, 1:end-1) + theta_u(:, 2:end)) / 2;
    T2 = T - 1;

    %% 3) Build distance matrix for linearly spaced intracranial sensors
    spacing = 3.5;   % mm

    idx = (1:N)';
    Dmat = abs(idx - idx.') * spacing;
    Dmat(1:N+1:end) = 0;

    %% 4) Parameterization
    % Unknowns:
    %   omega_i for i=1..N
    %   A0, T0, sigma, l
    %
    % x = [omega_1 ... omega_N, A0, T0, sigma, l]'

    % Initial guesses
    omega0 = mean(dTheta, 2);  % rough estimate from average derivative
    omega0 = min(max(omega0, omega_lo), omega_hi);

    A0_init = 1;
    T0_init = floor(T2/2);
    sigma_init = max(5, T2/10);
    l_init = median(Dmat(Dmat > 0));

    x0 = [omega0; A0_init; T0_init; sigma_init; l_init];

    % Bounds
    lb = [omega_lo * ones(N,1);  0;    1;      1e-3; 1e-3];
    ub = [omega_hi * ones(N,1);  Inf;  T2;     Inf;  Inf];

    %% 5) Define residual function
    % Residual stack over all channels and time points
    residual_fun = @(x) kuramoto_residual(x, dTheta, theta_mid, Dmat, fs);

    %% 6) Solve nonlinear least squares
    opts = optimoptions('lsqnonlin', ...
        'Display', 'iter', ...
        'MaxIterations', 200, ...
        'MaxFunctionEvaluations', 5e4, ...
        'FunctionTolerance', 1e-8, ...
        'StepTolerance', 1e-8);

    [xhat, resnorm, residual, exitflag, output] = lsqnonlin( ...
        residual_fun, x0, lb, ub, opts);

    %% 7) Unpack results
    omega = xhat(1:N);
    A0    = xhat(N+1);
    T0    = xhat(N+2);
    sigma = xhat(N+3);
    l     = xhat(N+4);

    params = struct();
    params.omega = omega;
    params.A0 = A0;
    params.T0 = T0;
    params.sigma = sigma;
    params.l = l;
    params.resnorm = resnorm;

    % Effective coupling matrix at the peak of A(t): A(T0)=A0
    G = exp(-Dmat ./ l);
    K_eff = A0 * G;
    K_eff(1:N+1:end) = 0;   % zero self-coupling

    fitinfo = struct();
    fitinfo.exitflag = exitflag;
    fitinfo.output = output;
    fitinfo.residual = residual;

end

function r = kuramoto_residual(x, dTheta, theta_mid, Dmat, fs)
% Residual vector for structured Kuramoto model.

    N = size(dTheta, 1);
    T2 = size(dTheta, 2);

    omega = x(1:N);
    A0    = x(N+1);
    T0    = x(N+2);
    sigma = x(N+3);
    l     = x(N+4);

    % Time index vector for midpoint samples
    t = (1:T2);

    % Gaussian envelope over time
    A_t = A0 * exp(-0.5 * ((t - T0) ./ sigma).^2);   % 1 x T2

    % Distance kernel
    G = exp(-Dmat ./ l);
    G(1:N+1:end) = 0;   % no self-coupling

    % Predicted derivatives
    dTheta_pred = zeros(N, T2);

    for i = 1:N
        % Construct spatially weighted sine sum
        S = zeros(1, T2);
        for j = 1:N
            if j == i
                continue;
            end
            S = S + G(i,j) * sin(theta_mid(j,:) - theta_mid(i,:));
        end

        dTheta_pred(i,:) = omega(i) + A_t .* S;
    end

    % Residuals stacked into a column vector
    r = (dTheta_pred - dTheta);
    r = r(:);

    % Optional ridge penalty on omega deviation from mid-range or on A0, etc.
    % You can uncomment and tune if desired.
    %
    % lambda_omega = 1e-4;
    % omega_ref = 2*pi*20;
    % r_reg = sqrt(lambda_omega) * (omega - omega_ref);
    % r = [r; r_reg];
end

function metrics = compute_fit_metrics(y_true, y_pred, num_params)
% y_true, y_pred: vectors of observed and predicted derivatives
% num_params: number of fitted parameters in the model

    y_true = y_true(:);
    y_pred = y_pred(:);

    residuals = y_true - y_pred;
    n = numel(y_true);

    RSS = sum(residuals.^2);
    TSS = sum((y_true - mean(y_true)).^2);

    metrics.RSS = RSS;
    metrics.MSE = RSS / n;
    metrics.RMSE = sqrt(metrics.MSE);
    metrics.R2 = 1 - RSS / TSS;

    % Information criteria
    metrics.AIC = 2 * num_params + n * log(RSS / n);
    metrics.BIC = num_params * log(n) + n * log(RSS / n);
end