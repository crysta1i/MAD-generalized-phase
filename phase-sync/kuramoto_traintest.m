function [K, omega, fitinfo] = kuramoto_traintest(allses_angle_cts, ch_names, fs)
% FIT_KURAMOTO_K
% Fit an unstructured Kuramoto model across multiple trials:
%
%   dtheta_i/dt = omega_i + sum_j K_ij * sin(theta_j - theta_i)
%
% Trial handling:
%   - Split trials into training and test sets
%   - Fit parameters jointly on all training trials
%   - Evaluate on held-out test trials
%
% Inputs:
%   allses_angle_cts : (n_timepoints x n_contacts x n_trials) phase data
%   ch_names         : cell array of channel names, length N
%   fs               : sampling frequency in Hz
%
% Outputs:
%   K       : N x N coupling matrix
%   omega   : N x 1 natural frequencies
%   fitinfo : struct with train/test metrics and AIC/BIC

    N = numel(ch_names);
    [~, n_contacts, n_trials] = size(allses_angle_cts);

    if n_contacts ~= N
        error('Number of contacts in allses_angle_cts does not match length(ch_names).');
    end

    % Train/test split
    rng(1);
    testFrac = 0.2;
    nTest = max(1, round(testFrac * n_trials));
    perm = randperm(n_trials);
    testTrials = perm(1:nTest);
    trainTrials = perm(nTest+1:end);

    fprintf('Unstructured model: %d training trials, %d test trials.\n', ...
        numel(trainTrials), numel(testTrials));

    % Frequency bounds
    f_lo = 5;
    f_hi = 40;
    omega_lo = 2*pi*f_lo;
    omega_hi = 2*pi*f_hi;

    omega = zeros(N, 1);
    K = zeros(N, N);

    trainResidualsAll = [];
    trainYAll = [];
    trainYhatAll = [];

    testResidualsAll = [];
    testYAll = [];
    testYhatAll = [];

    % Fit one regression per row i using stacked training trials
    for i = 1:N
        [Xtrain, ytrain] = build_unstructured_design(allses_angle_cts, trainTrials, fs, i);

        % Initial regularized fit with bounds on omega
        lb = [omega_lo; -Inf(N,1)];
        ub = [omega_hi;  Inf(N,1)];

        lambda = 0.1;
        reg_mask = [0; ones(N,1)];

        A_aug = [Xtrain; sqrt(lambda) * diag(reg_mask)];
        y_aug = [ytrain; zeros(N+1,1)];

        xhat = lsqlin(A_aug, y_aug, [], [], [], [], lb, ub, []);

        omega(i) = xhat(1);
        K(i,:)   = xhat(2:end)';

        % Training predictions for this row
        yhat_train = Xtrain * xhat;
        trainResidualsAll = [trainResidualsAll; (ytrain - yhat_train)]; %#ok<AGROW>
        trainYAll = [trainYAll; ytrain]; %#ok<AGROW>
        trainYhatAll = [trainYhatAll; yhat_train]; %#ok<AGROW>

        % Test predictions for this row
        [Xtest, ytest] = build_unstructured_design(allses_angle_cts, testTrials, fs, i);
        yhat_test = Xtest * xhat;
        testResidualsAll = [testResidualsAll; (ytest - yhat_test)]; %#ok<AGROW>
        testYAll = [testYAll; ytest]; %#ok<AGROW>
        testYhatAll = [testYhatAll; yhat_test]; %#ok<AGROW>
    end

    % Metrics
    trainRSS = sum(trainResidualsAll.^2);
    trainMSE = trainRSS / numel(trainResidualsAll);
    trainRMSE = sqrt(trainMSE);
    trainR2 = 1 - trainRSS / sum((trainYAll - mean(trainYAll)).^2);

    testRSS = sum(testResidualsAll.^2);
    testMSE = testRSS / numel(testResidualsAll);
    testRMSE = sqrt(testMSE);
    testR2 = 1 - testRSS / sum((testYAll - mean(testYAll)).^2);

    % Parameter count:
    % N omegas + N*N K entries
    % If you want to exclude self-coupling from parameter count, you can use N + N*(N-1)
    k = N + N*N;

    nTrainObs = numel(trainResidualsAll);
    AIC = 2*k + nTrainObs * log(trainRSS / nTrainObs);
    BIC = k * log(nTrainObs) + nTrainObs * log(trainRSS / nTrainObs);

    fitinfo = struct();
    fitinfo.trainTrials = trainTrials;
    fitinfo.testTrials = testTrials;
    fitinfo.trainRSS = trainRSS;
    fitinfo.trainMSE = trainMSE;
    fitinfo.trainRMSE = trainRMSE;
    fitinfo.trainR2 = trainR2;
    fitinfo.testRSS = testRSS;
    fitinfo.testMSE = testMSE;
    fitinfo.testRMSE = testRMSE;
    fitinfo.testR2 = testR2;
    fitinfo.AIC = AIC;
    fitinfo.BIC = BIC;
    fitinfo.k = k;
    fitinfo.trainResiduals = trainResidualsAll;
    fitinfo.testResiduals = testResidualsAll;
end


function [X, y] = build_unstructured_design(allses_angle_cts, trialIdx, fs, i)
% Build stacked design matrix and response vector for channel i across trials.

    [~, N, ~] = size(allses_angle_cts);

    Xcells = {};
    ycells = {};

    for r = 1:numel(trialIdx)
        tr = trialIdx(r);

        theta = squeeze(allses_angle_cts(:,:,tr))';  % N x T
        theta_u = unwrap(theta, [], 2);
        dTheta = diff(theta_u, 1, 2) * fs;           % N x (T-1)
        theta_mid = (theta_u(:,1:end-1) + theta_u(:,2:end)) / 2;

        T2 = size(dTheta, 2);

        y_r = dTheta(i,:)';                          % T2 x 1

        A = zeros(T2, N+1);
        A(:,1) = 1;                                  % intercept = omega_i
        for j = 1:N
            A(:,j+1) = sin(theta_mid(j,:) - theta_mid(i,:))';
        end

        Xcells{end+1} = A; %#ok<AGROW>
        ycells{end+1} = y_r; %#ok<AGROW>
    end

    X = vertcat(Xcells{:});
    y = vertcat(ycells{:});
end