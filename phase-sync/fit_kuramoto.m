% % EXAMPLE USAGE
% ch_names = cell(1, numel(cur_elec_contact_names));
% for i = 1:numel(cur_elec_contact_names)
%     ch_names{i} = cur_elec_contact_names(i);
% end
% trial    = 1;
% fs       = 2048;           

% K = fit_kuramoto(ch_names, trial, fs);

% % Visualise the coupling matrix
% figure;
% imagesc(K);
% colorbar;
% title('Kuramoto Coupling Matrix K_{ij}');
% xlabel('Source channel j');
% ylabel('Target channel i');
% xticklabels(ch_names);  yticklabels(ch_names);


function K = fit_kuramoto(ch_names, trial, fs)
% FIT_KURAMOTO  Estimate the N×N Kuramoto coupling matrix from EEG phase angles.
%
% Inputs:
%   ch_names : cell array of channel name strings, length N  e.g. {'Fp1','Fp2',...}
%   trial    : trial number passed through to get_phase_angles()
%   fs       : sampling frequency in Hz (used to convert diff() to rad/s)
%
% Output:
%   K        : N×N matrix of coupling constants K_ij  (rad/s)
%              K(i,j) is the influence of channel j on channel i.
%
% Model:
%   dtheta_i/dt = omega_i + sum_j  K_ij * sin(theta_j - theta_i)
%
% This is linearised in [omega_i, K_i1, ..., K_iN] and solved per row i
% via ordinary least squares.

    N = numel(ch_names);

    %% ── 1. Load phase-angle time series ──────────────────────────────────
    % theta: N × T  matrix  (radians)
    fprintf('Loading phase angles for %d channels, trial %d ...\n', N, trial);
    theta_cell = cell(N, 1);
    for i = 1 : N
        theta_cell{i} = get_phase_angles(ch_names{i}, i, trial);   % 1×T or T×1
    end

    % Unify orientation → N × T
    for i = 1 : N
        theta_cell{i} = theta_cell{i}(:)';          % force row vector
    end
    T_lengths = cellfun(@numel, theta_cell);
    if ~all(T_lengths == T_lengths(1))
        error('Channels have different numbers of time points.');
    end
    T = T_lengths(1);
    theta = vertcat(theta_cell{:});                 % N × T

    %% ── 2. Compute dtheta/dt via finite differences ──────────────────────

    % diff() gives T-1 differences; we use mid-points of theta for the RHS.
    % dTheta : N × (T-1),  units = radians/sample
    % Multiply by fs to get rad/s.
    dTheta = diff(theta, 1, 2) * fs;               % N × (T-1)

     % Handle phase-wrapping -- transitions from pi -> -pi
    idx_wrap = find(dTheta < -6);
    for idx = idx_wrap
        dTheta(idx) = (pi - theta(idx)) + (pi - abs(theta(idx+1)));
    end

    % Use the average of consecutive samples as the time point for the RHS
    theta_mid = (theta(:, 1:end-1) + theta(:, 2:end)) / 2;  % N × (T-1)
    T2 = T - 1;                                    % number of usable time points

    %% ── 3. Build the design matrix and solve one linear system per channel ─
    %
    % For channel i, the Kuramoto equation at each time t is:
    %
    %   dTheta_i(t) = omega_i * 1  +  K_i1*sin(theta_1(t)-theta_i(t))
    %                              +  K_i2*sin(theta_2(t)-theta_i(t))
    %                              +  ...
    %                              +  K_iN*sin(theta_N(t)-theta_i(t))
    %
    % Unknowns per row: x_i = [omega_i, K_i1, ..., K_iN]'   (N+1 unknowns)
    % Design matrix A_i : T2 × (N+1)
    %   col 1        → ones (for omega_i)
    %   col j+1      → sin(theta_j(t) - theta_i(t))   for j = 1..N
    %
    % Solve:  A_i * x_i  ≈  dTheta_i   (least squares)

    omega = zeros(N, 1);
    K     = zeros(N, N);

    for i = 1 : N
        % RHS: derivative of channel i  (T2 × 1)
        y = dTheta(i, :)';

        % Design matrix  (T2 × N+1)
        A = zeros(T2, N + 1);
        A(:, 1) = 1;                                % constant term → omega_i
        for j = 1 : N
            A(:, j + 1) = sin(theta_mid(j, :) - theta_mid(i, :))';
        end

        % Least-squares solution  (N+1 × 1)
        % Uses MATLAB's backslash which applies QR/SVD as appropriate
        x = A \ y;

        omega(i)   = x(1);
        K(i, :)    = x(2 : end)';                  % 1 × N  →  i-th row of K
    end

    %% ── 4. Report ─────────────────────────────────────────────────────────
    fprintf('\nEstimated natural frequencies omega (rad/s):\n');
    for i = 1 : N
        fprintf('  %s : %.4f\n', ch_names{i}, omega(i));
    end

    fprintf('\nCoupling matrix K (%dx%d) estimated successfully.\n', N, N);
end