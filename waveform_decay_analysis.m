% waveform_decay_analysis.m

% 1. for each alignment type, find time in trialavg of peak of interest (line 19 ref_peak_time) 
% 2. on single trials, find the peak closest to that time
% 3. check if at that time, peak has heavier right tail
%       - should this be amplitude based (how long it takes to get to same amplitude on left/right side) 
%         or phase based (how long it takes to get to same phase/trough on left/right side of the peak)?

% output proportion of trials on which this peak has a heavier right tail

addpath('/media/Data/Human_Intracranial_MAD/analysis/TravWaves/Code')
addpath('/media/Data/Human_Intracranial_MAD/_toolbox')

subject_ID = 'EMU038'; reference = 'neighbor_average'; cur_letter = "N"; 
phase_or_amp = "phase";
time_or_area = "area";
alignments = ["first_unique_attribute", "second_unique_attribute","third_unique_attribute","fourth_unique_attribute"];
[~, cur_elec_contact_names] = get_single_probe_contacts(reference, subject_ID, cur_letter);
ref_peak_time = 300; % will be different for different probes/contacts

% Using function #1
cnums = 10:-1:8;
name1 = cur_elec_contact_names(cnums(1)); name2 = cur_elec_contact_names(cnums(end));
fprintf("Proportion of trials where %s --> %s (ascending order of %s [%s])\n", name1, name2, time_or_area, phase_or_amp)
for alignment = alignments
    proportion_ordered_trials = check_contact_order(subject_ID, reference, cur_letter, cnums, alignment, ref_peak_time, phase_or_amp, time_or_area);
    fprintf("%s: %f \n", alignment, proportion_ordered_trials)
end

% % Using function #2
% contact = 10; ch_name = cur_elec_contact_names(contact);
% for alignment = alignments
%     [prop_heavy_right, diffs_right_left] = compute_prop_heavy_right(subject_ID, reference, cur_letter, contact, alignment, ref_peak_time, phase_or_amp, time_or_area);

%     fprintf("Proportion of %s trials in %s with heavier right tail (%s, %s): %f\n", alignment, ch_name, phase_or_amp, time_or_area, prop_heavy_right)
%     fprintf("Mean difference (right - left): %f\n", mean(diffs_right_left))
% end

% --------- Function definition -----------
function proportion_ordered_trials = check_contact_order(subject_ID, reference, cur_letter, cnums, alignment, ref_peak_time, phase_or_amp, time_or_area)
    % INPUTS:
    % - cnums: e.g., 7:9 or 9:-1:7 (input in the order in which the right AUC/decay time should be checked --
    % 7:9 means ct7 < ct8 < ct9 is considered "ordered" and 9:-1:7 means ct9 < ct8 < ct7 is desired order)
    
    % OUTPUT:
    % - proportion of trials where contacts' right AUC/decay time is ordered in accordance with contact position on probe
    
    all_ct_stat = []; 
    for contact = cnums
        [~, ~, AUC_time_right] = compute_prop_heavy_right(subject_ID, reference, cur_letter, contact, alignment, ref_peak_time, phase_or_amp, time_or_area);
        all_ct_stat = [all_ct_stat AUC_time_right];
    end
    signed_diffs = sign(diff(all_ct_stat, 1, 2));
    ordered_trials = zeros(height(all_ct_stat),1);
    for trial = 1:height(all_ct_stat)
        if all(signed_diffs(trial,:) == 1)
            ordered_trials(trial) = 1;
        end
    end
    proportion_ordered_trials = nnz(ordered_trials) / numel(ordered_trials);
end

function [proportion_heavy_right, diffs_right_left, AUC_time_right] = compute_prop_heavy_right(subject_ID, reference, cur_letter, contact, alignment, ref_peak_time, phase_or_amp, time_or_area)

% Inputs
% - contact: index in cur_elec_contact_ind corresponding to the contact 
%   on which this proportion shouuld be computed
% - ref_peak_time: time (in samples) of peak of interest in trial-avg
% - phase_or_amp: string, either "phase" or "amp", indicating whether 
%   to do phase-based or amplitude-based left/right tail comparison
% - time_or_area: string, either "time" or "area", indicating whether to 
%   define "heavy" as long decay time or large area under curve 

% Output:
% - proportion of single trials with a heavier right tail
% - diffs_right_left: list of time differences (in samples, not seconds) between right decay 
%   time and left rise time (N_trials x 1 array), right_diff - left_diff
% - AUC_time_right: list of right AUC (if "area") or right decay times (if "time")-- N_trials x 1 array

addpath('/media/Data/Human_Intracranial_MAD/analysis/TravWaves/Code')
addpath('/media/Data/Human_Intracranial_MAD/analysis/TravWaves/Code/ERP')
addpath('/media/Data/Human_Intracranial_MAD/analysis/TravWaves/Code/generalized-phase/analysis')
[~, data_base_dir, ~, ~, num_sessions, Fs, ~] = tw_setup(subject_ID, reference);
[cur_elec_contact_ind, ~] = get_single_probe_contacts(reference, subject_ID, cur_letter);
ch_num = cur_elec_contact_ind(contact);

if strcmp(alignment, "trialstart") || strcmp(alignment, "anticipation")
    win = [Fs/2 0];
else
    win = [0 Fs/2];
end

inspect_alignments = ["first_unique_attribute", "second_unique_attribute", "third_unique_attribute", "fourth_unique_attribute", "inspection", "single_opt_first_inspection", "full_single_opt_info"];
if ismember(alignment, inspect_alignments)
    align_type = "inspection";
end

if ~strcmpi(time_or_area, "area") && ~strcmpi(time_or_area, "time")
    error("invalid argument for time_or_area")
elseif ~strcmpi(phase_or_amp, "amp") && ~strcmpi(phase_or_amp, "phase")
    error("invalid argument for phase_or_amp")
elseif ref_peak_time > abs(win(2) - win(1))
    error("ref_peak_time out of bounds")
end

% % Compute ref_peak_time: time of highest peak in trial_avg
% [xgp_avg_allct1, avg_events_allct1, ~, ~, ~, ~, ~, ~] = plot_save_trial_avg(subject_ID, reference, align_type, alignment, alignment, cur_letter, false, false, false);
% trialavg = avg_events_allct1(:, contact);
% trialavg_angle = angle(xgp_avg_allct1(:, contact));
% [~, time_highest_peak] = max(trialavg);

% if strcmpi(phase_or_amp, "phase") % use generalized phase
%     peak_times = find(diff(sign(trialavg_angle)) > 0); % find times when this time series crosses 0 phase, take the time closest to time_highest_peak
%     [~, idx_ref] = min(abs(peak_times - time_highest_peak));
%     ref_peak_time = peak_times(idx_ref);
% elseif strcmpi(phase_or_amp, "amp") % use local min/max
%     ref_peak_time = time_highest_peak;
% else
%     error("cannot handle argument for phase_or_amp param")
% end

count = 0; num_events = 0; diffs_right_left = []; AUC_time_right = [];
right_empty = 0; left_empty = 0;
for sesnum = 1:num_sessions
    if strcmp(reference,'Ground')
        load(sprintf('%s/%s/%s_MAD_SES%d_Setup.mat',data_base_dir,subject_ID,subject_ID,sesnum), "filters", "trial_times", "trial_words");
    else
        load(sprintf('%s/%s/%s_MAD_SES%d_Setup_%s.mat',data_base_dir,subject_ID,subject_ID,sesnum,reference), "filters", "trial_times", "trial_words");
    end

    [align_times,~] = get_align_times(filters, trial_times, trial_words, alignment);
    align_times(isnan(align_times)) = [];
    align_times = round(align_times*Fs);  

    if strcmp(subject_ID,'EMU001') && strcmp(align_type, 'inspection')
        [align_times_inspect,~] = get_align_times(filters, trial_times, trial_words, '2opt_2att_inspection');
        align_times_inspect(isnan(align_times_inspect)) = [];
        align_times_inspect = round(align_times_inspect*Fs); 
        align_times = intersect(align_times, align_times_inspect);
    end

    if strcmp(reference,'Ground')
        load(sprintf('%s/%s/separate_channel_files/%s_MAD_SES%d_ch%03d.mat',data_base_dir,subject_ID,subject_ID,sesnum,ch_num),"data"); 
    elseif strcmp(reference,'neighbor_average')            
        load(sprintf('%s/%s/separate_channel_files/%s/%s_MAD_SES%d_ch%03d.mat',data_base_dir,subject_ID,reference,subject_ID,sesnum,ch_num),"data");
    end

    [b,a] = butter( 4, [5 40] ./ (Fs/2) ); filtered = filtfilt( b, a, data ); data = filtered;

    for event = 1:numel(align_times)
        % local max method to identify peaks
        data_event_padded = data(round(align_times(event)-win(1)-Fs):round(align_times(event)+win(2)-1)+Fs);
        data_event = data(round(align_times(event)-win(1)):round(align_times(event)+win(2)-1));
        zscored = (data_event - mean(data_event)) / std(data_event);

        if strcmpi(phase_or_amp, "amp")

            peak_times_event = find(diff(sign(diff(data_event))) < 0);
            trough_times_event = find(diff(sign(diff(data_event))) > 0);
            [~, peak_idx] = min(abs(peak_times_event - ref_peak_time));
            peak_time = peak_times_event(peak_idx);
            
            trough_before_time = max(trough_times_event(trough_times_event < peak_time));
            trough_after_time = min(trough_times_event(trough_times_event > peak_time));

            % % ------- amplitude -------- % OPTION 1
            % ref_trough_amp = max(zscored(trough_before_time), zscored(trough_after_time));
            % times_ref_amp = find(diff(sign(zscored - ref_trough_amp)) ~= 0);
            %     % find time at which data crosses from being less than the ref_trough_amp to being greater
            % time_ref_right = min(times_ref_amp(times_ref_amp > peak_time));
            % time_ref_left = max(times_ref_amp(times_ref_amp < peak_time));
            
            % -------- local-min trough -------- % OPTION 2
            time_ref_left = trough_before_time; 
            time_ref_right = trough_after_time;

            if strcmpi(time_or_area, "area")
                % % ------- fixed time -------- % OPTION 3
                % time_ref_right = peak_time + min(trough_after_time - peak_time, peak_time - trough_before_time);
                % time_ref_left = peak_time - min(trough_after_time - peak_time, peak_time - trough_before_time);
                % % ---------------------------

                left = sum(zscored((time_ref_left:peak_time-1)));
                right = sum(zscored(peak_time+1:time_ref_right));
            elseif strcmpi(time_or_area, "time")
                left = peak_time - time_ref_left;
                right = time_ref_right - peak_time;
            end

        elseif strcmpi(phase_or_amp, "phase")

            xgp_event_padded = generalized_phase_vector(data_event_padded, Fs, 5);
            angle_event = angle(xgp_event_padded(1+Fs:1+Fs+abs(win(2)-win(1))-1));

            peak_times_event = find(diff(sign(angle_event)) > 0);
            trough_times_event = find(diff(sign(angle_event)) < 0);
            [~, peak_idx] = min(abs(peak_times_event - ref_peak_time));
            peak_time = peak_times_event(peak_idx);

            trough_before_time = max([max([win(1)+1; peak_time - 100]); trough_times_event(trough_times_event < peak_time)]); %max([win(1)+1; peak_time - 100]); 
            trough_after_time = min([min([peak_time + 100; win(2)-1]); trough_times_event(trough_times_event > peak_time)]); %; min([peak_time + 100; win(2)-1]);

            if strcmpi(time_or_area, "area")
                left = sum(zscored((trough_before_time:peak_time-1)));
                right = sum(zscored(peak_time+1:trough_after_time));
            elseif strcmpi(time_or_area, "time")
                left = peak_time - trough_before_time;
                right = trough_after_time - peak_time;
            end

        end

        if right > left
            count = count + 1;
        end 

        if isempty(right) 
            right_empty = right_empty + 1;
        elseif isempty(left)
            left_empty = left_empty + 1;
        end

        num_events = num_events + 1;
        diffs_right_left = [diffs_right_left; right-left];
        AUC_time_right = [AUC_time_right; right];
    end
end
proportion_heavy_right = count / num_events;

end