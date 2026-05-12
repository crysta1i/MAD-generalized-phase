function [s_plus, s_minus, s_recon] = tw_fourier_reconstruct(a_plus, a_minus, k_vals, v_est, x, t)
% tw_fourier_reconstruct (real)
%   Traveling wave fourier basis RECONSTRUCTION
%   Reconstruct real-valued right- and left-going traveling-wave components.
%
% INPUTS
%   a_plus, a_minus : complex vectors
%       Coefficients from tw_fourier_decomp.
%   k_vals : vector
%       Spatial wavenumbers corresponding to coefficients.
%   v_est : float
%       Estimated velocity of the traveling wave.
%   x : vector (N_x)
%       Spatial positions at which to reconstruct.
%   t : vector (N_t)
%       Time samples at which to reconstruct.
%
% OUTPUTS
%   s_plus : [N_x x N_t] array
%       Real-valued right-going reconstructed component.
%   s_minus : [N_x x N_t] array
%       Real-valued left-going reconstructed component.
%   s_recon : [N_x x N_t] array
%       Total reconstruction (s_plus + s_minus).

[X, T] = ndgrid(x, t);

s_plus  = zeros(size(X));
s_minus = zeros(size(X));

for ik = 1:length(k_vals)
    k = k_vals(ik);
    % Right-going basis: cos(kx - v k t) and sin(kx - v k t)
    phase_plus  = k*X - v_est*k*T;
    s_plus = s_plus + real(a_plus(ik)) * cos(phase_plus) ...
                    - imag(a_plus(ik)) * sin(phase_plus);

    % Left-going basis: cos(kx + v k t) and sin(kx + v k t)
    phase_minus = k*X + v_est*k*T;
    s_minus = s_minus + real(a_minus(ik)) * cos(phase_minus) ...
                      - imag(a_minus(ik)) * sin(phase_minus);
end

s_recon = s_plus + s_minus;

end
