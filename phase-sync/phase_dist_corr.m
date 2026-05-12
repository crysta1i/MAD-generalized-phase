function [corr_mat, pv_mat] = phase_dist_corr(ch_nums, subwin_st, subwin_end, allses_angle_cts)

% PHASE-DIST CORRELATION
% Compute circular-linear correlation at each timept (in specified subwindow) 
% between phase angle and distance/position of contacts on a single probe

% INPUTS
% - cur_letter: string, probe name
% - chnums: array of integers (indices in cur_elec_contact_ind)
% - alignments: string array of specific alignments (e.g., ["single_opt_first_inspection", "full_single_opt_info"])
% - subwin_st: number of samples (after behavioral alignment) to start computing 
% - allses_angle_cts: output from function compute_angle_ts(subject_ID, reference, cur_letter, ch_nums, alignment, subwin_st, subwin_end)
% OUTPUTS
% - corr_mat: each column is time series of correlation coefficients for one trial; column number = trial number
% - pv_mat: p-values (size of matrix corresponds to that of corr_mat)

addpath('/media/Data/Human_Intracranial_MAD/analysis/TravWaves/Code')
addpath('/media/Data/Human_Intracranial_MAD/analysis/PAC_code/CF_Coupling/generalized phase/circstat-matlab')

corr_mat = zeros(subwin_end - subwin_st + 1,0); % timepoints x trials (eventually avearge rows to avg across trials)
pv_mat = zeros(subwin_end - subwin_st + 1,0); 
% report number of significant trials? but need to do mc correction

% iterate through trials (3rd dimension of angles_cts)
for trial = 1:size(allses_angle_cts, 3)
    corr_vec = zeros(subwin_end - subwin_st + 1, 1);
    pv_vec = zeros(subwin_end - subwin_st + 1, 1);

    for time = subwin_st:subwin_end
        [rho, pval] = circ_corrcl(allses_angle_cts(time - subwin_st + 1,:,trial), 1:numel(ch_nums));
        corr_vec(time - subwin_st + 1) = rho;
        pv_vec(time - subwin_st + 1) = pval;
    end
    corr_mat(:,end+1) = corr_vec;
    pv_mat(:,end+1) = pv_vec;
end

end