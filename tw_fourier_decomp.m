function [v_est, a_plus, a_minus, k_vals, omega_vals] = tw_fourier_decomp(signal_xt, dx, dt, doWindow)

% tw_fourier_decomp
%   Decompose a spatio-temporal signal into traveling-wave components.
%   traveling wave fourier basis DECOMPOSITION (ChatGPT)
%
% INPUTS
%   signal_xt : [N_x x N_t] array
%       Signal measured at N_x spatial positions and N_t time points.
%   dx : float
%       Spatial spacing between sensors.
%   dt : float
%       Temporal spacing between samples.
%   doWindow : logical
%       If true, apply Hann window in both dimensions.
%
% OUTPUTS
%   v_est : float
%       Estimated wave speed (absolute value).
%   a_plus, a_minus : complex vectors
%       Coefficients along the right-going (+) and left-going (-) ridges.
%   k_vals : vector
%       Spatial wavenumbers (rad/unit length).
%   omega_vals : vector
%       Temporal frequencies (rad/unit time).

[N_x, N_t] = size(signal_xt);
L = N_x * dx;
T = N_t * dt;

sig = signal_xt;

% Windowing
if doWindow
    win_x = hann(N_x);
    win_t = hann(N_t).';
    sig = sig .* (win_x * win_t);
end

% 2D FFT with fftshift
S = fftshift(fft2(sig));

k_vals = 2*pi*fftshift(fftfreq(N_x, dx));
omega_vals = 2*pi*fftshift(fftfreq(N_t, dt));

% Build grids
[K, W] = meshgrid(k_vals, omega_vals);
K = K.';  % match MATLAB indexing (rows = x, cols = t)
W = W.';

% Candidate velocities
v_candidates = linspace(0, 0.9*dx/dt, 200);
ridge_energy = zeros(size(v_candidates));

tol = 2*pi/T;  % frequency bin tolerance

for iv = 1:length(v_candidates)
    v = v_candidates(iv);
    mask_plus  = abs(W -  v*K) < tol;
    mask_minus = abs(W +  v*K) < tol;
    ridge_energy(iv) = sum(abs(S(mask_plus)).^2 + abs(S(mask_minus)).^2);
end

[~, imax] = max(ridge_energy);
v_est = v_candidates(imax);

% Extract coefficients along ridges
a_plus  = zeros(size(k_vals));
a_minus = zeros(size(k_vals));
for ik = 1:length(k_vals)
    k = k_vals(ik);
    [~, idx_plus]  = min(abs(omega_vals - v_est*k));
    [~, idx_minus] = min(abs(omega_vals + v_est*k));
    a_plus(ik)  = S(ik, idx_plus);
    a_minus(ik) = S(ik, idx_minus);
end

end

% --- Helper function: fftfreq equivalent ---
function f = fftfreq(N, d)
    % Frequencies in cycles per unit (like numpy.fft.fftfreq)
    if mod(N,2)==0
        f = [0:(N/2-1) -N/2:-1] / (N*d);
    else
        f = [0:((N-1)/2) -((N-1)/2):-1] / (N*d);
    end
end
