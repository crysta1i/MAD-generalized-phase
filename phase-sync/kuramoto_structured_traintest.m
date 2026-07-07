function [params, omega, K_eff, fitinfo] = kuramoto_structured_traintest(allses_angle_cts, ch_names, fs)
% KURAMOTO_STRUCTURED
% Fit a Kuramoto model with structured, time-varying coupling across multiple trials:
%
%   dtheta_i/dt = omega_i + A(t) * sum_j exp(-d_ij/l) * sin(theta_j - theta_i)
%
% where
%   A(t) = A0 * exp(-(t - T0)^2 / (2*sigma^2))
%
% Trial handling:
%   - Split trials into training and test sets
%   - Fit parameters jointly on all training trials
%   - Evaluate prediction error on held-out test trials
%
% Inputs:
%   allses_angle_cts : (n_timepoints x n_contacts x n_trials) phase data
%   ch_names         : cell array of channel name strings, length N
%   fs               : sampling frequency in Hz
%
% Outputs:
%   params           : struct with fitted parameters
%   omega            : N x 1 fitted natural frequencies (rad/s)
%   K_eff            : N x N effective coupling matrix at peak A(t)
%   fitinfo          : struct with optimization info and prediction errors

    %% Settings
    N = numel(ch_names);
    [n_timepoints, n_contacts, n_trials] = size(allses_angle_cts);

    if n_contacts ~= N
        error('Number of contacts in allses_angle_cts does not match length(ch_names).');
    end

    % Train/test split
    rng(1); % reproducible split
    testFrac = 0.2;
    nTest = max(1, round(testFrac * n_trials));
    perm = randperm(n_trials);
    testTrials = perm(1:nTest);
    trainTrials = perm(nTest+1:end);

    fprintf('Structured model: %d training trials, %d test trials.\n', ...
        numel(trainTrials), numel(testTrials));

    %% Frequency bounds for omega in rad/s
    f_lo = 5;
    f_hi = 40;
    omega_lo = 2*pi*f_lo;
    omega_hi = 2*pi*f_hi;

    %% Build distance matrix for linearly spaced contacts
    spacing = 3.5e-3; % meters
    idx = (1:N)';
    Dmat = abs(idx - idx.') * spacing;
    Dmat(1:N+1:end) = 0;

    %% Build training data stacked across trials
    trainData = build_trial_stack(allses_angle_cts, trainTrials, fs);

    % trainData fields:
    %   dTheta_all   : N x Tstack
    %   theta_mid_all : N x Tstack
    %   trial_time    : 1 x Tstack, local time index within each trial (1..T2)
    %   trial_id      : 1 x Tstack, identifies which trial each sample came from

    dTheta_train = trainData.dTheta_all;
    theta_train   = trainData.theta_mid_all;
    tIdx_train    = trainData.trial_time;

    T2 = size(dTheta_train, 2);

    %% Initial guesses
    omega0 = mean(dTheta_train, 2);
    omega0 = min(max(omega0, omega_lo), omega_hi);

    A0_init = 1;
    T0_init = round(T2 / 2);
    sigma_init = max(5, T2 / 10);
    l_init = 5 * spacing;  % a few contacts
    x0 = [omega0; A0_init; T0_init; sigma_init; l_init];

    %% Bounds
    lb = [omega_lo * ones(N,1);  0;     1;      1e-3; 1e-6];
    ub = [omega_hi * ones(N,1);  Inf;   T2;     Inf;   Inf];

    %% Nonlinear least squares fit on training trials
    residual_fun = @(x) structured_residual(x, dTheta_train, theta_train, Dmat, tIdx_train);

    opts = optimoptions('lsqnonlin', ...
        'Display', 'iter', ...
        'MaxIterations', 200, ...
        'MaxFunctionEvaluations', 5e4, ...
        'FunctionTolerance', 1e-8, ...
        'StepTolerance', 1e-8);

    [xhat, resnorm, residual, exitflag, output] = lsqnonlin( ...
        residual_fun, x0, lb, ub, opts);

    %% Unpack fitted parameters
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
    params.trainTrials = trainTrials;
    params.testTrials = testTrials;

    %% Effective coupling at peak A(t)
    G = exp(-Dmat ./ l);
    G(1:N+1:end) = 0;
    K_eff = A0 * G;

    %% Evaluate prediction error on train and test sets
    trainMetrics = evaluate_structured_fit(allses_angle_cts, trainTrials, fs, omega, A0, T0, sigma, l, Dmat);
    testMetrics  = evaluate_structured_fit(allses_angle_cts, testTrials,  fs, omega, A0, T0, sigma, l, Dmat);

    %% AIC / BIC on training data
    % Number of parameters:
    % N omegas + 4 global parameters
    k = N + 4;

    % Use training residual sum of squares for IC
    RSS_train = sum(trainMetrics.residuals.^2);
    nTrainObs = numel(trainMetrics.residuals);

    AIC = 2*k + nTrainObs * log(RSS_train / nTrainObs);
    BIC = k * log(nTrainObs) + nTrainObs * log(RSS_train / nTrainObs);

    %% Package fit info
    fitinfo = struct();
    fitinfo.exitflag = exitflag;
    fitinfo.output = output;
    fitinfo.residual = residual;
    fitinfo.trainMetrics = trainMetrics;
    fitinfo.testMetrics = testMetrics;
    fitinfo.AIC = AIC;
    fitinfo.BIC = BIC;
    fitinfo.k = k;

end


function data = build_trial_stack(allses_angle_cts, trialIdx, fs)
% Build stacked derivative and midpoint-phase arrays over a set of trials.
%
% Output:
%   dTheta_all   : N x Tstack
%   theta_mid_all : N x Tstack
%   trial_time   : 1 x Tstack local time index (1..T2)
%   trial_id     : 1 x Tstack trial label

    [n_timepoints, N, ~] = size(allses_angle_cts);

    dTheta_cells = {};
    thetaMid_cells = {};
    trialTime_cells = {};
    trialId_cells = {};

    for r = 1:numel(trialIdx)
        tr = trialIdx(r);

        theta = squeeze(allses_angle_cts(:,:,tr))';  % N x T
        theta_u = unwrap(theta, [], 2);

        dTheta = diff(theta_u, 1, 2) * fs;           % N x (T-1)
        theta_mid = (theta_u(:,1:end-1) + theta_u(:,2:end)) / 2;

        T2 = size(dTheta, 2);

        dTheta_cells{end+1} = dTheta; %#ok<AGROW>
        thetaMid_cells{end+1} = theta_mid; %#ok<AGROW>
        trialTime_cells{end+1} = 1:T2; %#ok<AGROW>
        trialId_cells{end+1} = tr * ones(1, T2); %#ok<AGROW>
    end

    data.dTheta_all = cat(2, dTheta_cells{:});
    data.theta_mid_all = cat(2, thetaMid_cells{:});
    data.trial_time = cat(2, trialTime_cells{:});
    data.trial_id = cat(2, trialId_cells{:});
end


function r = structured_residual(x, dTheta, theta_mid, Dmat, tIdx)
% Residual vector for structured Kuramoto model across stacked trials.

    N = size(dTheta, 1);
    T2 = size(dTheta, 2);

    omega = x(1:N);
    A0    = x(N+1);
    T0    = x(N+2);
    sigma = x(N+3);
    l     = x(N+4);

    G = exp(-Dmat ./ l);
    G(1:N+1:end) = 0;

    A_t = A0 * exp(-0.5 * ((tIdx - T0) ./ sigma).^2);  % 1 x T2

    dTheta_pred = zeros(N, T2);

    for i = 1:N
        S = zeros(1, T2);
        for j = 1:N
            if j == i
                continue;
            end
            S = S + G(i,j) * sin(theta_mid(j,:) - theta_mid(i,:));
        end
        dTheta_pred(i,:) = omega(i) + A_t .* S;
    end

    r = (dTheta_pred - dTheta);
    r = r(:);
end


function metrics = evaluate_structured_fit(allses_angle_cts, trialIdx, fs, omega, A0, T0, sigma, l, Dmat)
% Evaluate structured model on selected trials.

    [~, N, ~] = size(allses_angle_cts);

    residuals_all = [];
    y_all = [];
    yhat_all = [];

    G = exp(-Dmat ./ l);
    G(1:N+1:end) = 0;

    for r = 1:numel(trialIdx)
        tr = trialIdx(r);

        theta = squeeze(allses_angle_cts(:,:,tr))';  % N x T
        theta_u = unwrap(theta, [], 2);
        dTheta = diff(theta_u, 1, 2) * fs;           % N x (T-1)
        theta_mid = (theta_u(:,1:end-1) + theta_u(:,2:end)) / 2;

        T2 = size(dTheta, 2);
        tIdx = 1:T2;
        A_t = A0 * exp(-0.5 * ((tIdx - T0) ./ sigma).^2);

        dTheta_pred = zeros(N, T2);

        for i = 1:N
            S = zeros(1, T2);
            for j = 1:N
                if j == i
                    continue;
                end
                S = S + G(i,j) * sin(theta_mid(j,:) - theta_mid(i,:));
            end
            dTheta_pred(i,:) = omega(i) + A_t .* S;
        end

        residuals = dTheta_pred - dTheta;

        residuals_all = [residuals_all; residuals(:)]; %#ok<AGROW>
        y_all = [y_all; dTheta(:)]; %#ok<AGROW>
        yhat_all = [yhat_all; dTheta_pred(:)]; %#ok<AGROW>
    end

    RSS = sum(residuals_all.^2);
    TSS = sum((y_all - mean(y_all)).^2);

    metrics.RSS = RSS;
    metrics.MSE = RSS / numel(residuals_all);
    metrics.RMSE = sqrt(metrics.MSE);
    metrics.R2 = 1 - RSS / TSS;
    metrics.residuals = residuals_all;
    metrics.y = y_all;
    metrics.yhat = yhat_all;
    metrics.nObs = numel(residuals_all);
end