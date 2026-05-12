function allses_angle_cts = compute_angle_ts(subject_ID, reference, cur_letter, ch_nums, alignment, subwin_st, subwin_end)

% ANGLE_TS
% Compute circular-linear correlation at each timept (in specified subwindow) 
% between phase angle and distance/position of contacts on a single probe

% INPUTS
% - cur_letter: string, probe name
% - chnums: array of integers (indices in cur_elec_contact_ind)
% - alignments: string array of specific alignments (e.g., ["single_opt_first_inspection", "full_single_opt_info"])
% - subwin_st: number of samples (after behavioral alignment) to start computing 
% OUTPUTS
% - allses_angle_cts: 3D array -- time series of angles for all contacts in ch_nums and all trials (time x cnums x trials)

addpath('/media/Data/Human_Intracranial_MAD/analysis/TravWaves/Code')
addpath('/media/Data/Human_Intracranial_MAD/_toolbox')

% TODO: function to compute xgp on single trials

[~, data_base_dir, ~, ~, num_sessions, Fs, ~] = tw_setup(subject_ID, reference);
[cur_elec_contact_ind, ~] = get_single_probe_contacts(reference, subject_ID, cur_letter);
win = [0 Fs/2];

allses_angle_cts = zeros(subwin_end - subwin_st + 1, numel(ch_nums), 0);
for sesnum = 1:numel(num_sessions)
    if strcmp(reference,'Ground')
        load(sprintf('%s/%s/%s_MAD_SES%d_Setup.mat',data_base_dir,subject_ID,subject_ID,sesnum), "filters", "trial_times", "trial_words");
    else
        load(sprintf('%s/%s/%s_MAD_SES%d_Setup_%s.mat',data_base_dir,subject_ID,subject_ID,sesnum,reference), "filters", "trial_times", "trial_words");
    end

    [align_times,~] = get_align_times(filters, trial_times, trial_words, alignment);
    align_times(isnan(align_times)) = [];
    align_times = round(align_times*Fs);  

    if strcmp(subject_ID,'EMU001') && strcmp(alignment, 'inspection')
        [align_times_inspect,~] = get_align_times(filters, trial_times, trial_words, '2opt_2att_inspection');

        align_times_inspect(isnan(align_times_inspect)) = [];
        align_times_inspect = round(align_times_inspect*Fs); 

        align_times = intersect(align_times, align_times_inspect);
    end

    angle_cts = zeros(subwin_end - subwin_st + 1, numel(ch_nums), numel(align_times)); % 3D: timepts x channels x trials
    for ct = 1:numel(ch_nums)
        cnum = ch_nums(ct);
        contact = cur_elec_contact_ind(cnum);
        if strcmp(reference,'Ground')
            load(sprintf('%s/%s/separate_channel_files/%s_MAD_SES%d_ch%03d.mat',data_base_dir,subject_ID,subject_ID,sesnum,contact),"data"); 
        elseif strcmp(reference,'neighbor_average')            
            load(sprintf('%s/%s/separate_channel_files/%s/%s_MAD_SES%d_ch%03d.mat',data_base_dir,subject_ID,reference,subject_ID,sesnum,contact),"data");
        end
            
        [b,a] = butter( 4, [5 40] ./ (Fs/2) ); filtered = filtfilt( b, a, data ); data = filtered;

        % for every event, clip the signal (with padding) and add to data_cts
        %time_by_trials = zeros(subwin_end - subwin_st + 1, numel(ch_nums));
        for event = 1:numel(align_times)
            data_padded = data(round(align_times(event)-win(1)) - Fs:round(align_times(event)+win(2)-1) + Fs);
            xgp = generalized_phase_vector(data_padded, Fs, 5);
            xgp = xgp(1+Fs:1+Fs+win(2)-1);
            angle_vec = angle(xgp);

            angle_cts(:, ct, event) = angle_vec(subwin_st:subwin_end);
        end
    end
    allses_angle_cts = cat(3, allses_angle_cts, angle_cts);

end % end loop through sessions

end