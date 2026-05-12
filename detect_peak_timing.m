% detect_peak_timing.m

% Find patterns in timing of peak/trough of spatially neighboring contacts on an electrode

% Save results in folder: phase_progression_detection
% Create a cell array for each electrode:
%      -Column 1: Event number (inspection #)
%      -Column 2: Subwindow start (subwin_start in samples; e.g., 150 = 150th time sample in the event window)
%      -Column 8: Peaks or troughs detected to progress? (0 if peaks, 3.1415 if troughs) 
%      -Column 9: Speed of wave propagation (units = contacts / millisecond)
%      -Column 10: Increasing contacts [ct1, ct2] -- timings/latencies increase starting from ct1 and ending at ct2
%                 (empty array if no increasing subsquence >= 4)
%      -Column 11: Decreasing contacts [ct1, ct2] -- latencies decrease starting from ct1 and ending at ct2 
%                 (empty array if no decreasing subsequence >= 4)
%      -Column 12: 100x2 logical vector indicating circular correlation coeffcient and 
%       significance (1=yes, 0=no) of it for each timepoint in the subwindow
%      -Column 13: Direction of timing progression (increasing = 1, decreasing = -1, both = 0)
%      -Column 14: Corr_sig -- 100x2 double table, column 1 = circular correlation coefficients 
%                 at each time point in sfprintfubwindow, column 2 = significance (0/1)
%      *NOTE: contact_num = 1 is actually the largest number contact (most lateral)
% Convert cell array to table before saving (cell2table)

%% 1. Using fixed subwindow sizes (100 samples) and overlap lengths (50 samples)

subject_IDs = {'EMU024','EMU025','EMU030','EMU038','EMU039','EMU041'};
reference = 'Ground'; sesnum = 1; Fs = 2048;
alignment = 'inspection'; % nonmotor % inspection % trialstart
win = [0 1024]; %[0 1024] [1024 0]
datapath = '/media/Data/Human_Intracranial_MAD/';
data_base_dir = sprintf('%s1_formatted',datapath);


for subject_num = numel(subject_IDs) % loop through subjects
    subject_ID = subject_IDs{subject_num};

    % all setup files
    final_out_dir = sprintf('%s/generalized_phase/waveform_plots/trial_averaged/%s',out_dir, alignment);
    if ~exist(final_out_dir, 'dir'), mkdir(final_out_dir); end 
    
    table_out_dir = sprintf('%s/phase_progression_detection',final_out_dir);
    if ~exist(table_out_dir, 'dir'), mkdir(table_out_dir); end

    load(sprintf('%s/%s/%s_MAD_SES%d_Setup.mat',data_base_dir,subject_ID,subject_ID,sesnum));
    %if strcmp(subject_ID, 'EMU025'), elec_name = elec_name(1:206); end
    
    % load brainstorm file and outer contacts
    cTable = readtable(sprintf('/media/Data/Human_Intracranial_MAD/analysis/Location_Consistency/3D Coordinates/%s.tsv',subject_ID), "FileType","delimitedtext",'Delimiter', '\t');
    channel_names_bs = cTable.Channel; % extract "Channels" variable in table

    % get electrode letters
    elec_letters = [];
    for i = 1:numel(channel_names_bs)
        cur_name = channel_names_bs{i}; 
        if contains(cur_name, "'")
            if ~ismember(cur_name(1:2),elec_letters)
                elec_letters = [elec_letters; convertCharsToStrings(cur_name(1:2))];
            end
        else 
            if ~ismember(cur_name(1),elec_letters)
                elec_letters = [elec_letters; convertCharsToStrings(cur_name(1))];
            end
        end
        
    end

    % get align times
    [align_times,trial_numbers] = get_align_times(filters, trial_times, trial_words, alignment);
    remove_ind = isnan(align_times);
    align_times(isnan(align_times)) = []; trial_numbers(remove_ind) = [];
    align_times = round(align_times*Fs);

    if numel(align_times) < 30
        max_events = numel(align_times); 
    else
        max_events = 30;
    end

    % TODO: include more columns in phase_prog_table marking the type ofinspection:
        % 6/15: win/loss, amt/prob, first att, novel/repeat (binary 0/1 coding)

    % 'win_domain','loss_domain','amount','probability','novel_attribute_inspection','repeated_inspection',...
    % 'first_unique_attribute','second_unique_attribute','third_unique_attribute','fourth_unique_attribute'

    if strcmp(alignment,'inspection')
        win_align_times = get_align_times(filters, trial_times, trial_words, 'win_domain');
        win_align_times(isnan(win_align_times)) = [];
        win_align_times = round(win_align_times*Fs);
            
        %loss_align_times = get_align_times(filters, trial_times, trial_words, 'loss_domain');
        amt_align_times = get_align_times(filters, trial_times, trial_words, 'amount');
        amt_align_times(isnan(amt_align_times)) = [];
        amt_align_times = round(amt_align_times*Fs);

        %prob_align_times = get_align_times(filters, trial_times, trial_words, 'probability');
        first_unique_align_times = get_align_times(filters, trial_times, trial_words, 'first_unique_attribute');
        first_unique_align_times(isnan(first_unique_align_times)) = [];
        first_unique_align_times = round(first_unique_align_times*Fs);

        fourth_unique_align_times = get_align_times(filters, trial_times, trial_words, 'fourth_unique_attribute');
        fourth_unique_align_times(isnan(fourth_unique_align_times)) = [];
        fourth_unique_align_times = round(fourth_unique_align_times*Fs);
        
        novel_align_times = get_align_times(filters, trial_times, trial_words, 'novel_attribute_inspection');
        novel_align_times(isnan(novel_align_times)) = [];
        novel_align_times = round(novel_align_times*Fs);
        %repeat_align_times = get_align_times(filters, trial_times, trial_words, 'repeated_inspection');
    end

    for elec_num = 1:numel(elec_letters) % loop through electrodes
        
        cur_letter = elec_letters(elec_num);

        phase_prog_table = cell(0,14);
        for event = 1:numel(align_times) % loop through events

            win_bool = ismember(align_times(event), win_align_times);
            amt_bool = ismember(align_times(event), amt_align_times);
            first_bool = ismember(align_times(event), first_unique_align_times);
            fourth_bool = ismember(align_times(event), fourth_unique_align_times);
            novel_bool = ismember(align_times(event), novel_align_times); 

            xgp_cur_elec = []; % each column = 1 contact on current electrode, each row = 1 timept
            j = 1;
            if ~contains(cur_letter, "'")
                while ~contains(elec_name{j},cur_letter) || contains(elec_name{j},"'"), j = j+1; end
            else
                while ~contains(elec_name{j},cur_letter), j = j+1; end
            end
            cur_name = elec_name{j};
            cur_elec_contact_names = []; % stores name of contacts on cur electrode
            cur_elec_contact_ind = []; % stores indices of contacts on cur elec
            while strcmp(cur_name(1:2),cur_letter) || strcmp(cur_name(1),cur_letter)
                cur_elec_contact_names = [cur_elec_contact_names; convertCharsToStrings(cur_name)];
                cur_elec_contact_ind = [cur_elec_contact_ind; j];
                load(sprintf('%s/extracted_phase/xgp_%s_ses%02d_ch%03d.mat',out_dir,cur_name,sesnum,j));
                
                xgp_cur_contact = xgp((round(align_times(event)-win(1)):round(align_times(event)+win(2)-1)));
                xgp_cur_elec = [xgp_cur_elec xgp_cur_contact];
                
                j = j+1; 
                if j > numel(elec_name), break; end
                cur_name = elec_name{j};
            end
            
            for subwin_start = [1, 50:50:900] % loop through subwindows of 100 samples within event windows
                subwindow = [subwin_start subwin_start + 100];
                angles_cur_elec = angle(xgp_cur_elec);
                angles_cur_elec = angles_cur_elec(subwindow(1):subwindow(2),:);
                for phase_type = [0 1] % check both peak 0 and trough pi at each subwindow
                    ref_phase = pi*phase_type;
                    peak_times = zeros(1,width(angles_cur_elec));
                    for e = 1:width(angles_cur_elec)
                        [peak1, peak_time] = min(abs(angles_cur_elec(:,e) - ref_phase));
                        peak_times(e) = peak_time;
                    end
                    % detect monotone sequences in this list
                    max_incr_count = 1; max_decr_count = 1;
                    max_incr_end = 1; max_decr_end = 1;
                    increase_count = 1; decrease_count = 1;
                    prev_latency = peak_times(1);
                    prev_diff = abs(peak_times(1) - peak_times(2));
                    idx = 2;
                    while idx <= numel(peak_times)
                        while idx <= numel(peak_times) && peak_times(idx) > prev_latency && abs(abs(peak_times(idx) - prev_latency) - prev_diff) < 15
                            if ~(peak_times(idx) < 100 && prev_latency > 1)
                                prev_latency = peak_times(idx); idx = idx + 1;
                                break;
                            end
                            increase_count = increase_count + 1;
                            prev_diff = abs(peak_times(idx) - prev_latency);
                            prev_latency = peak_times(idx);
                            idx = idx + 1;
                        end
                        if increase_count > max_incr_count
                            max_incr_count = increase_count; 
                            max_incr_end = idx;
                        end
                        increase_count = 1; % reset
                        prev_diff = abs(peak_times(idx) - peak_times(idx + 1));

                        while idx <= numel(peak_times) && peak_times(idx) < prev_latency && abs(abs(prev_latency - peak_times(idx)) - prev_diff) < 15
                            if ~(peak_times(idx) > 1 && prev_latency < 100)
                                prev_latency = peak_times(idx); idx = idx + 1;
                                break;
                            end
                            decrease_count = decrease_count + 1;
                            prev_diff = abs(peak_times(idx) - prev_latency);
                            prev_latency = peak_times(idx);
                            idx = idx + 1;
                        end
                        if decrease_count > max_decr_count
                            max_decr_count = decrease_count; 
                            max_decr_end = idx; % max_decr_start = max_decr_end - max_decr_count
                        end
                        decrease_count = 1; % reset
                        prev_diff = abs(peak_times(idx) - peak_times(idx + 1));

                        while idx <= numel(peak_times) && peak_times(idx) == prev_latency
                            prev_latency = peak_times(idx);
                            idx = idx + 1;
                        end
 
                    end
                    if max_incr_count >= 4
                        incr_speed = (max_incr_count - 1) / ((peak_times(max_incr_end - 1) - peak_times(max_incr_end - max_incr_count)) * (1000/2048));
                    end
                    if max_decr_count >= 4
                        decr_speed = (max_decr_count - 1) / ((peak_times(max_decr_end - 1) - peak_times(max_decr_end - max_decr_count)) * (1000/2048));
                    end
                    % Ensure that peak timings follow roughly linear progression:
                    incr_cond = max_incr_count >= 4; %&& range(diff(peak_times((max_incr_end - max_incr_count):(max_incr_end - 1)))) <= 20;
                    decr_cond = max_decr_count >= 4; %&& range(diff(peak_times((max_decr_end - max_decr_count):(max_decr_end - 1)))) <= 20;
                    if incr_cond || decr_cond 
                        if incr_cond && decr_cond
                            direction = 0; % both increasing and decreasing latencies
                            ct_incr = [max_incr_end - max_incr_count, max_incr_end - 1];
                            ct_decr = [max_decr_end - max_decr_count, max_decr_end - 1];
                            if max_incr_count >= max_decr_count
                                speed = incr_speed;
                            else
                                speed = decr_speed;
                            end
                        elseif incr_cond
                            ct_incr = [max_incr_end - max_incr_count, max_incr_end - 1];
                            ct_decr = [];
                            direction = 1; % increasing latencies
                            speed = incr_speed;
                        else
                            ct_decr = [max_decr_end - max_decr_count, max_decr_end - 1];
                            ct_incr = [];
                            direction = -1; % decreasing latencies
                            speed = decr_speed;
                        end
                        % create corr_sig
                        corr_sig = zeros(100,2); % zeros(33,2) TODO
                        for kk = 1:100 % length of subwindow % 3:3:100
                            pl = angles_cur_elec(kk,:); % j is time point index
                            [ cc, pv ] = circ_corrcl(pl, 1:numel(cur_elec_contact_names)); % all in a row so exact coordinates don't matter
                            corr_sig(kk,:) = [cc, pv < 0.05];
                        end
                        % TODO: approximate linearity check before including in table
                        phase_prog_table(end + 1, :) = {event, subwin_start, win_bool, amt_bool, first_bool, fourth_bool, novel_bool, ref_phase, speed, ct_incr, ct_decr, peak_times, direction, corr_sig};
                    end
                end % end loop through phase type
            end % end loop through subwindows
        end % end loop through events
        
        phase_prog_table = cell2table(phase_prog_table,"VariableNames",["Event","Subwin St","Win","Amt","1st","4th","Novel","Peak/Trough","Speed (ct/ms)","Incr Contacts","Decr Contacts","Latencies","Direction","Corr"]);

        fname_table = sprintf('%s/%s_%s_phase_prog_table.mat', table_out_dir, cur_letter, alignment); 
        save(fname_table, "phase_prog_table")
        disp("saved")
    
    end % end loop through electrodes
end % end loop through subjects

%%
column5 = phase_prog_table.Corr;
majority_sig_count = 0;
for ii=1:numel(column5)
    cur_corr_sig = column5{ii}; 
    if nnz(cur_corr_sig(:,2)) > 50
        majority_sig_count = majority_sig_count + 1; 
    end 
end

%% Post-processing:
% Filtering rows of the table

idx_filt = phase_prog_table{:,"Subwin St"} >= 1 & phase_prog_table{:,"Subwin St"} <= 200; %& phase_prog_table{:,"Direction"} > -1;

    % includes time filter:
idx_filt = phase_prog_table{:,"Subwin St"} == 250; %& phase_prog_table{:,"Amt"} == 0 & phase_prog_table{:,"Novel"} == 1; % & phase_prog_table{:,"1st"} == 0
filtTable = phase_prog_table(idx_filt,:);

    % doesn't include time filter:
idx_event_filt = phase_prog_table{:,"Amt"} == 0 & phase_prog_table{:,"1st"} == 0 & phase_prog_table{:,"Novel"} == 1;
filtEventTable = phase_prog_table(idx_event_filt,:);

load("/media/Data/Human_Intracranial_MAD/analysis/Cross_Freq_Coupling/EMU024/Ground/synchrony_data/generalized_phase/phase_progression_detection/F_inspection_phase_prog_table.mat")
load("/media/Data/Human_Intracranial_MAD/analysis/Cross_Freq_Coupling/EMU024/Ground/synchrony_data/generalized_phase/phase_progression_detection/F_trialstart_phase_prog_table.mat")

%% count number of rows, num incr and num decr

subject_IDs = {'EMU024','EMU025','EMU030','EMU038','EMU039','EMU041'};
for subject_num = 1:numel(subject_IDs) % for each subject
    subject_ID = subject_IDs{subject_num};

    % get list of electrode letters
    cTable = readtable(sprintf('/media/Data/Human_Intracranial_MAD/analysis/Location_Consistency/3D Coordinates/%s.tsv',subject_ID), "FileType","delimitedtext",'Delimiter', '\t');
    channel_names_bs = cTable.Channel; 
    elec_letters = [];
    for i = 1:numel(channel_names_bs)
        cur_name = channel_names_bs{i};
        if contains(cur_name, "'")
            if ~ismember(cur_name(1:2),elec_letters)
                elec_letters = [elec_letters; convertCharsToStrings(cur_name(1:2))];
            end
        else
            if ~ismember(cur_name(1),elec_letters)
                elec_letters = [elec_letters; convertCharsToStrings(cur_name(1))];
            end
        end
    end
    
    fprintf("\n-------- %s --------", subject_ID)
    for el = 1:numel(elec_letters) % for each electrode letter
        cur_letter = elec_letters(el);
        
        % load inspection pahse prog table
        load(sprintf("/media/Data/Human_Intracranial_MAD/analysis/Cross_Freq_Coupling/%s/Ground/synchrony_data/generalized_phase/phase_progression_detection/%s_inspection_phase_prog_table.mat",subject_ID, cur_letter))
        filter_early_idx = phase_prog_table{:,"Subwin St"} < 500;
        filter_early = phase_prog_table(filter_early_idx,:);
        count_early = height(filter_early);
        count_early_incr = nnz(filter_early.Direction == 1);
        count_early_decr = nnz(filter_early.Direction == -1);

        filter_late_idx = phase_prog_table{:,"Subwin St"} >= 500;
        filter_late = phase_prog_table(filter_late_idx,:);
        count_late = height(filter_late);
        count_late_incr = nnz(filter_late.Direction == 1);
        count_late_decr = nnz(filter_late.Direction == -1);
        
        % count_inspect = height(phase_prog_table);
        % count_inspect_incr = nnz(phase_prog_table.Direction == 1);
        % count_inspect_decr = nnz(phase_prog_table.Direction == -1);

        % load trialstart phase prog table
        % load(sprintf("/media/Data/Human_Intracranial_MAD/analysis/Cross_Freq_Coupling/%s/Ground/synchrony_data/generalized_phase/phase_progression_detection/%s_trialstart_phase_prog_table.mat",subject_ID, cur_letter))
        % count_trialst = height(phase_prog_table);
        % count_trialst_incr = nnz(phase_prog_table.Direction == 1);
        % count_trialst_decr = nnz(phase_prog_table.Direction == -1);

        fprintf("\n%s early: %d detections (%d incr, %d dec)", cur_letter, count_early, count_early_incr, count_early_decr)
        fprintf("\n%s late: %d detections (%d incr, %d dec)", cur_letter, count_late, count_late_incr, count_late_decr)
        fprintf("\n")

    end
end

%% 2. Using variable subwindow sizes based on peak/trough times of reference contact(s)
% Reference1 = contact 1 (innermost contact)
% Reference2 = contact numel(cur_elec_contact_names) (outermost contact)

% 6/16/25 NEW IDEA
% Use all of the contacts as reference: find the time of the peaks for each
% of the contacts on an electrode -> find where most of these peaks are
% clustered in time -> find the *first* and *last* peaks in each cluster to
% define subwindows 
% (how to determine whether a peak belongs to "this" cluster or "next"
% cluster of peaks?)
% Next, detect progressions of peaks in time, and detect progressions of
% other per-peak phases in time (e.g., phase = 0.1, -0.1, etc). Check if
% the progression of other phases follows a *similar* pattern as the
% progression of peaks in time

% Also, I want to test whether 2+ signals have similar *shape* (same
% waveform but simply translated across time/space) -- I can use a
% similarity measure like xcorr? or check similarity of peak progressions

subject_IDs = {'EMU024','EMU025','EMU030','EMU038','EMU039'};
reference = 'Ground'; sesnum = 1; Fs = 2048;
alignment = 'inspection'; % nonmotor
win = [0 1024]; %[0 1024] [1024 0]
datapath = '/media/Data/Human_Intracranial_MAD/';
data_base_dir = sprintf('%s1_formatted',datapath);

for subject_num = 1:numel(subject_IDs) % loop through subjects
    subject_ID = subject_IDs{subject_num};

    % all setup files
    final_out_dir = sprintf('/media/Data/Human_Intracranial_MAD/analysis/Cross_Freq_Coupling/%s/%s/synchrony_data/generalized_phase',subject_ID, reference);
    table_out_dir = sprintf('%s/phase_progression_detection',final_out_dir);
    if ~exist(table_out_dir, 'dir'), mkdir(table_out_dir); end

    load(sprintf('%s/%s/%s_MAD_SES%d_Setup.mat',data_base_dir,subject_ID,subject_ID,sesnum));
    %if strcmp(subject_ID, 'EMU025'), elec_name = elec_name(1:206); end
    
    % load brainstorm file and outer contacts
    cTable = readtable(sprintf('/media/Data/Human_Intracranial_MAD/analysis/Location_Consistency/3D Coordinates/%s.tsv',subject_ID), "FileType","delimitedtext",'Delimiter', '\t');
    channel_names_bs = cTable.Channel; % extract "Channels" variable in table

    % get electrode letters
    elec_letters = [];
    for i = 1:numel(channel_names_bs)
        cur_name = channel_names_bs{i};
        if contains(cur_name, "'")
            if ~ismember(cur_name(1:2),elec_letters)
                elec_letters = [elec_letters; convertCharsToStrings(cur_name(1:2))];
            end
        else 
            if ~ismember(cur_name(1),elec_letters)
                elec_letters = [elec_letters; convertCharsToStrings(cur_name(1))];
            end
        end
        
    end

    % get align times
    [align_times,trial_numbers] = get_align_times(filters, trial_times, trial_words, alignment);
    remove_ind = isnan(align_times);
    align_times(isnan(align_times)) = []; trial_numbers(remove_ind) = [];
    align_times = round(align_times*Fs);

    for elec_num = 1:numel(elec_letters) % loop through electrodes
        
        cur_letter = elec_letters(elec_num);

        phase_prog_table = cell(0,7);
        for event = 1:30 %numel(align_times) % loop through events
            
            xgp_cur_elec = []; % each column = 1 contact on current electrode, each row = 1 timept
            j = 1;
            if ~contains(cur_letter, "'")
                while ~contains(elec_name{j},cur_letter) || contains(elec_name{j},"'"), j = j+1; end
            else
                while ~contains(elec_name{j},cur_letter), j = j+1; end
            end
            cur_name = elec_name{j};
            cur_elec_contact_names = []; % stores name of contacts on cur electrode
            cur_elec_contact_ind = []; % stores indices of contacts on cur elec
            while strcmp(cur_name(1:2),cur_letter) || strcmp(cur_name(1),cur_letter)
                cur_elec_contact_names = [cur_elec_contact_names; convertCharsToStrings(cur_name)];
                cur_elec_contact_ind = [cur_elec_contact_ind; j];
                load(sprintf('%s/extracted_phase/xgp_%s_ses%02d_ch%03d.mat',final_out_dir,cur_name,sesnum,j));
                
                xgp_cur_contact = xgp((round(align_times(event)-win(1)):round(align_times(event)+win(2)-1)));
                xgp_cur_elec = [xgp_cur_elec xgp_cur_contact];
                
                j = j+1; 
                if j > numel(elec_name), break; end
                cur_name = elec_name{j};
            end
            angles_cur_elec = angle(xgp_cur_elec); % t rows, c columns (c = # of contacts)
            
            % Make list of subwindow starts by finding peaks/troughs
                % To find peaks, find 0 crossings in angles_cur_elec(:,4) and angles_cur_elec(:,end - 3)
            idx_cross_peak_first = diff(angles_cur_elec(:,4) > 0) ~= 0;
            idx_cross_peak_last = diff(angles_cur_elec(:,end - 3) > 0) ~= 0;
                % To find troughs, same process as finding peaks but first subtract pi from all angles?
                % but need phase angles to wrap around!
            rotated_angles_first = wrapToPi(angles_cur_elec(:,4) - pi);
            rotated_angles_last = wrapToPi(angles_cur_elec(:, end-3) - pi);
            idx_cross_trough_first = diff(rotated_angles_first > 0) ~= 0;
            idx_cross_trough_last = diff(rotated_angles_last > 0) ~= 0;

            % note: 1st subwindow should span from beginning of window to
            % midway between 1st and 2nd peaks; 2nd subwindow spans from
            % midway between 1st and 2nd peaks to midway btwn 2nd and 3rd

            % Find midpoint between each neighboring pair of entries in idx_cross_peak_first
            subwin_st_peaks_first = movmean(idx_cross_peak_first, 2);
            subwin_st_troughs_first = movmean(idx_cross_trough_first, 2); % first = first half, last = last half
            subwin_st_peaks_last = movmean(idx_cross_peak_last,2);
            subwin_st_troughs_last = movmean(idx_cross_trough_last,2);
                % TODO: is it important to track which contact was used as the reference?
                % TODO: should we make these windows overlap? add one additional window in
                % between each pair of these?
            subwin_st_troughs = [1; subwin_st_troughs_first; subwin_st_troughs_last];
            subwin_st_peaks = [1; subwin_st_peaks_first; subwin_st_peaks_last];

            for phase_type = [0 1] % check both peak 0 and trough pi at each subwindow
                ref_phase = pi*phase_type;
                if phase_type == 0
                    subwin_starts = subwin_st_peaks;
                else
                    subwin_starts = subwin_st_troughs;
                end
                
                for subwin_start = subwin_starts
                    subwindow = [subwin_start subwin_start + 100];
                    angles_cur_elec = angle(xgp_cur_elec);
                    angles_cur_elec = angles_cur_elec(subwindow(1):subwindow(2),:);
                    
                    peak_times = zeros(1,width(angles_cur_elec));
                    for e = 1:width(angles_cur_elec)
                        [peak1, peak_time] = min(abs(angles_cur_elec(:,e) - ref_phase));
                        peak_times(e) = peak_time;
                    end
                    % detect monotone sequences in this list
                    max_incr_count = 1; max_decr_count = 1;
                    max_incr_end = 1; max_decr_end = 1;
                    increase_count = 1; decrease_count = 1;
                    prev_latency = peak_times(1);
                    idx = 2;
                    while idx <= numel(peak_times)
                        while idx <= numel(peak_times) && peak_times(idx) > prev_latency && peak_times(idx) < 100 && prev_latency > 1
                            increase_count = increase_count + 1;
                            prev_latency = peak_times(idx);
                            idx = idx + 1;
                        end
                        if increase_count > max_incr_count
                            max_incr_count = increase_count; 
                            max_incr_end = idx;
                        end
                        increase_count = 1; % reset

                        while idx <= numel(peak_times) && peak_times(idx) < prev_latency && peak_times(idx) > 1 && prev_latency < 100
                            decrease_count = decrease_count + 1;
                            prev_latency = peak_times(idx);
                            idx = idx + 1;
                        end
                        if decrease_count > max_decr_count
                            max_decr_count = decrease_count; 
                            max_decr_end = idx; % max_decr_start = max_decr_end - max_decr_count
                        end
                        decrease_count = 1; % reset

                        while idx <= numel(peak_times) && peak_times(idx) == prev_latency
                            prev_latency = peak_times(idx);
                            idx = idx + 1;
                        end
 
                    end
                    if max_incr_count >= 4 || max_decr_count >= 4
                        incr_speed = (max_incr_count - 1) / ((peak_times(max_incr_end - 1) - peak_times(max_incr_end - max_incr_count)) * (1000/2048));
                        decr_speed = (max_decr_count - 1) / ((peak_times(max_decr_end - 1) - peak_times(max_decr_end - max_decr_count)) * (1000/2048));
                        
                        if max_incr_count >=4 && max_decr_count >= 4
                            direction = 0; % both increasing and decreasing latencies
                            if max_incr_count >= max_decr_count
                                speed = incr_speed;
                            else
                                speed = decr_speed;
                            end
                        elseif max_incr_count >=4
                            direction = 1; % increasing latencies
                            speed = incr_speed;
                        else
                            direction = -1; % decreasing latencies
                            speed = decr_speed;
                        end
                        % create corr_sig
                        corr_sig = zeros(100,2); % zeros(33,2) TODO
                        for kk = 1:100 % length of subwindow % 3:3:100
                            pl = angles_cur_elec(kk,:); % j is time point index
                            [ cc, pv ] = circ_corrcl(pl, 1:numel(cur_elec_contact_names)); % all in a row so exact coordinates don't matter
                            corr_sig(kk,:) = [cc, pv < 0.05];
                        end

                        phase_prog_table(end + 1, :) = {event, subwin_start, direction, ref_phase, corr_sig, speed, peak_times};
                    end
                end % end loop through phase type
            end % end loop through subwindows
        end % end loop through events
        
        phase_prog_table = cell2table(phase_prog_table,"VariableNames",["Event","Subwindow Start","Direction","Peak/Trough","Corr","Speed (ct/ms)","Latencies"]);
        fname_table = sprintf('%s/%s_%s_refcontact_phaseprog_table.mat', table_out_dir, cur_letter, alignment); 
        save(fname_table, "phase_prog_table")
    
    end % end loop through electrodes
end % end loop through subjects

%% TEMP copy of section 2 of detect_waves_trialavg
reference = 'Ground';
alignments = ["inspection","single_opt_first_inspection","full_single_opt_info"];
datapath = '/media/Data/Human_Intracranial_MAD/';
data_base_dir = sprintf('%s1_formatted',datapath);
addpath('/media/Data/Human_Intracranial_MAD/analysis/PAC_code/CF_Coupling/generalized phase/circstat-matlab')

subject_IDs = {'EMU001','EMU024','EMU025','EMU030','EMU037','EMU038','EMU039','EMU041','EMU047','EMU051'};

for subject_num = 1 % [8 9 5] %1:numel(subject_IDs)
    subject_ID = subject_IDs{subject_num}; 
    
    if strcmp(subject_ID,'EMU001')
        Fs = 1000;
    else
        Fs = 2048;
    end
    
    out_dir = sprintf('/media/Data/Human_Intracranial_MAD/analysis/Cross_Freq_Coupling/%s/%s/synchrony_data/generalized_phase',subject_ID, reference);    

    switch subject_ID
        case 'EMU001'
            num_sessions = 3;
        case 'EMU024'
            num_sessions = 3;
        case 'EMU025'
            num_sessions = 2;
        case 'EMU030'
            num_sessions = 2;
        case 'EMU037'
            num_sessions = 4;
        case 'EMU038'
            num_sessions = 1;
        case 'EMU039'
            num_sessions = 4;
        case 'EMU041'
            num_sessions = 9;
        case 'EMU047'
            num_sessions = 1;
        case 'EMU051'
            num_sessions = 1;
    end

    load(sprintf('%s/%s/%s_MAD_SES%d_Setup.mat',data_base_dir,subject_ID,subject_ID,1));

    channel_names_bs = elec_name;

    % get elec letters
    elec_letters = [];
    for i = 1:numel(channel_names_bs)
        cur_name = channel_names_bs{i};
        if contains(cur_name, "'")
            if ~ismember(cur_name(1:2),elec_letters)
                elec_letters = [elec_letters; convertCharsToStrings(cur_name(1:2))];
            end
        else
            if ~ismember(cur_name(1),elec_letters)
                elec_letters = [elec_letters; convertCharsToStrings(cur_name(1))];
            end
        end
    end   
        
        for e = 1 %:numel(elec_letters) % TEMP 11/1/25
            cur_letter = elec_letters(e);
            load(sprintf('%s/fdr_ERP_sig_table_%s.mat', out_dir, cur_letter))

            j = 1;
            if ~contains(cur_letter, "'")
                while ~contains(elec_name{j},cur_letter) || contains(elec_name{j},"'"), j = j+1; end
            else
                while ~contains(elec_name{j},cur_letter), j = j+1; end
            end
            cur_name = elec_name{j};
            cur_elec_contact_names = []; % stores name of contacts on cur electrode
            cur_elec_contact_ind = []; % stores indices of contacts on cur elec
            
            while strcmp(cur_name(1:2),cur_letter) || strcmp(cur_name(1),cur_letter)
                cur_elec_contact_names = [cur_elec_contact_names; convertCharsToStrings(cur_name)];
                cur_elec_contact_ind = [cur_elec_contact_ind; j];
                
                j = j+1; 
                if j > numel(elec_name), break; end
                cur_name = elec_name{j};
            end

            ERP_peaks = ERP_table(ERP_table{:,"Peak/Trough"} == "peak",:);
            ERP_troughs = ERP_table(ERP_table{:,"Peak/Trough"} == "trough",:);
            % ERP_peaks = sortrows(ERP_peaks,4,"ascend"); % sort Idx column
            % ERP_troughs = sortrows(ERP_troughs,4,"ascend"); % sort Idx column

            % Create new table for saving traveling waves detected
            traveling_ERP_table = cell(0,10);
            

            for alignment = alignments

                per_contact_wave_counts = zeros(numel(cur_elec_contact_ind),1);
                per_contact_ERP_counts = zeros(numel(cur_elec_contact_ind),1);

                travel_ERP_count = 0;
                %fprintf("--------Proportion of ERPs in traveling waves--------\n")
                fprintf("%s %s: \n ", alignment, cur_letter)

                out_dir = sprintf('/media/Data/Human_Intracranial_MAD/analysis/Cross_Freq_Coupling/%s/%s/synchrony_data/generalized_phase',subject_ID, reference);  
                final_out_dir = sprintf('%s/waveform_plots/trial_averaged/%s',out_dir, alignment);
                data_out_dir = sprintf('%s/data_filtered_trial_avg/%s',out_dir, alignment);
            
                ERP_peaks_align = ERP_peaks(strcmp(ERP_peaks{:,"Alignment"},alignment),:);
                ERP_troughs_align = ERP_troughs(strcmp(ERP_troughs{:,"Alignment"},alignment),:);

                % Delete all rows in ERP_peaks_align and ERP_troughs_align
                % that are not significant
                ERP_peaks_align(ERP_peaks_align{:,"Sig"}==0,:) = [];
                ERP_troughs_align(ERP_troughs_align{:,"Sig"}==0,:) = [];

                % 10/30/25 - Create temporary new columns in ERP_peaks_align and ERP_troughs_align
                % (boolean) true = this row's ERP participates in traveling wave; 
                %           false = does not participate in traveling wave
                trav_pks_bool = zeros(height(ERP_peaks_align),1);
                trav_tghs_bool = zeros(height(ERP_troughs_align),1);

                %if isempty(ERP_peaks_align) && isempty(ERP_troughs_align), continue; end

                    if ~isempty(ERP_peaks_align)
                        ERP_peaks_align.IdxBefore = zeros(height(ERP_peaks_align),1);
                        ERP_peaks_align.IdxAfter = zeros(height(ERP_peaks_align),1);
                        for row = 1:height(ERP_peaks_align)
                            peak_idx = ERP_peaks_align{row,"Idx"};
                            name = ERP_peaks_align{row,"Contact"};
                            ct = find(strcmp(elec_name,name));
                            cnum = find(strcmp(cur_elec_contact_names,name));
                            per_contact_ERP_counts(cnum) = per_contact_ERP_counts(cnum) + 1;
                            
                            fname_data_avg = sprintf("%s/realavg_%s_%s_ch%03d.mat", data_out_dir, alignment, name, ct);
                            load(fname_data_avg, "avg_events")

                            fname_xgp_avg = sprintf("%s/centered_%s_xgp_%s_ch%03d.mat", data_out_dir, alignment, name, ct);
                            load(fname_xgp_avg, "xgp_avg")

                            % find all troughs (angle xgp method)
                            angle_posneg = xgp_avg;
                            angle_posneg(angle_posneg >= 0) = 1;
                            angle_posneg(angle_posneg < 0) = -1;
                            troughs_idx_align = find(diff(angle_posneg) < 0);
                            troughs_align = avg_events(troughs_idx_align);
                            
                            % %find all troughs (derivative method)
                            % first_deriv_align = diff(avg_events);
                            % first_deriv_align(first_deriv_align < 0) = -1;
                            % first_deriv_align(first_deriv_align > 0) = 1;
                            % troughs_idx_align = find(diff(first_deriv_align) > 0);
                            % troughs_align = avg_events(troughs_idx_align);
                            
                            % troughs before/after
                            idx_trough_before = troughs_idx_align(troughs_idx_align < peak_idx);
                            if ~isempty(idx_trough_before)
                                idx_trough_before = idx_trough_before(end);
                            else
                                idx_trough_before = 1;
                            end
                            ERP_peaks_align{row,"IdxBefore"} = idx_trough_before;
        
                            idx_trough_after = troughs_idx_align(troughs_idx_align > peak_idx);
                            if ~isempty(idx_trough_after)
                                idx_trough_after = idx_trough_after(1);
                            else
                                idx_trough_after = numel(avg_events);
                            end
                            ERP_peaks_align{row,"IdxAfter"} = idx_trough_after;
    
                            % Test each row: take [IdxBefore, IdxAfter] as the subwindow
                            subwin_start = idx_trough_before;
                            subwin_end = idx_trough_after;


                            xgp_avg_cur_elec = []; % each column = 1 contact on current electrode, each row = 1 timept
                            avg_cur_elec = [];
                            j = 1;
                            if ~contains(cur_letter, "'")
                                while ~contains(elec_name{j},cur_letter) || contains(elec_name{j},"'"), j = j+1; end
                            else
                                while ~contains(elec_name{j},cur_letter), j = j+1; end
                            end
                            cur_name = elec_name{j};
                            cur_elec_contact_names = []; % stores name of contacts on cur electrode
                            cur_elec_contact_ind = []; % stores indices of contacts on cur elec

                            while strcmp(cur_name(1:2),cur_letter) || strcmp(cur_name(1),cur_letter)
                                cur_elec_contact_names = [cur_elec_contact_names; convertCharsToStrings(cur_name)];
                                cur_elec_contact_ind = [cur_elec_contact_ind; j];

                                data_dir = sprintf('%s/data_filtered_trial_avg/%s',out_dir, alignment);
                                fname_xgp = sprintf("%s/centered_%s_xgp_%s_ch%03d.mat", data_dir, alignment, cur_name, j);
                                load(fname_xgp)
                                xgp_avg_cur_elec = [xgp_avg_cur_elec xgp_avg];

                                fname_data_avg = sprintf("%s/realavg_%s_%s_ch%03d.mat", data_dir, alignment, cur_name, j);
                                load(fname_data_avg) % loads avg_events
                                %avg_cur_elec = [avg_cur_elec avg_events];

                                j = j+1; 
                                if j > numel(elec_name), break; end
                                cur_name = elec_name{j};
                            end

                            angles_diff_sign_elec = angle(xgp_avg_cur_elec); % t rows, c columns (c = # of contacts)

                           % AT THIS BREAK POINT: check if angles_cur_elec has correct dimensions
                            angles_subwin = angles_diff_sign_elec(subwin_start:subwin_end,:);
                            %data_subwin = avg_cur_elec(subwin_starts(j):subwin_ends(j),:);
                    
                            % FIND TIME of PEAKS
                            ref_phase = 0;
                            peak_times = zeros(1,numel(cur_elec_contact_names)); % PEAKS
                            for k = 1:numel(cur_elec_contact_names)
                                [peak1, peak_time] = min(abs(angles_subwin(:,k) - ref_phase));
                                
                                % data_deriv = diff(data_subwin(:,k)); 
                                % data_posneg = data_deriv;
                                % data_posneg(data_posneg >= 0) = 1; data_posneg(data_posneg < 0) = -1;
                                % allpeaks_idx = find(diff(data_posneg) < 0) + 1; allpeaks = data_subwin(allpeaks_idx);
                                % [max_peak, max_time] = max(allpeaks);
                                % peak_time = allpeaks_idx(max_time);
                    
                                peak_times(k) = peak_time;
                            end
                            
                            % FIND MONOTONE SEQUENCES
                            %incr_wave_idx = []; decr_wave_idx = [];

                            max_incr_count = 1; max_decr_count = 1;
                            max_incr_end = 1; max_decr_end = 1; %max_incr_st = 1; max_decr_st = 1;
                            increase_count = 1; decrease_count = 1;
                            prev_latency = peak_times(1);
                            idx = 2;
                            while idx <= numel(peak_times)
                                while idx <= numel(peak_times) && peak_times(idx) > prev_latency 
                                    if ~(peak_times(idx) < subwin_end && prev_latency > 1)
                                        prev_latency = peak_times(idx); idx = idx + 1;
                                        break;
                                    end
                                    increase_count = increase_count + 1;
                                    prev_latency = peak_times(idx);
                                    idx = idx + 1;
                                end
                                if increase_count >= 4
                                    %incr_wave_idx = [incr_wave_idx; [idx-increase_count idx-1]];
                                    p = polyfit(1:increase_count,peak_times(idx-increase_count:idx-1),1);
                                    incr_speed = 3.5 / (p(1) * (1000/2048));
                                    circlin_corr = zeros(height(angles_subwin),2);
                                    for t = 1:height(angles_subwin)
                                        [cc, pv] = circ_corrcl(1:increase_count,angles_subwin(t,idx-increase_count:idx-1));
                                        circlin_corr(t,1) = cc; circlin_corr(t,2) = pv;
                                    end
                                    traveling_ERP_table(end + 1, :) = {alignment, "peak", subwin_start, subwin_end, [idx-increase_count idx-1], 1, peak_times, incr_speed, circlin_corr, ERP_table{row,1}};
                                    if (cnum >= idx-increase_count) & (cnum <= idx-1) % TEMP added 10/10/25
                                        travel_ERP_count = travel_ERP_count + 1;
                                        per_contact_wave_counts(cnum) = per_contact_wave_counts(cnum) + 1;
                                        trav_pks_bool(row) = 1;
                                    end
                                end
                                % if increase_count > max_incr_count
                                %     max_incr_count = increase_count; 
                                %     max_incr_end = idx;
                                % end
                                increase_count = 1; % reset
                               
                                while idx <= numel(peak_times) && peak_times(idx) < prev_latency 
                                    if ~(peak_times(idx) > 1 && prev_latency < subwin_end)
                                        prev_latency = peak_times(idx); idx = idx + 1;
                                        break;
                                    end
                                    decrease_count = decrease_count + 1;
                                    prev_latency = peak_times(idx);
                                    idx = idx + 1;
                                end

                                if decrease_count >= 4
                                    %decr_wave_idx = [decr_wave_idx; [idx-decrease_count idx-1]];
                                    p = polyfit(1:decrease_count,peak_times(idx-decrease_count:idx-1),1);
                                    decr_speed = 3.5 / (p(1) * (1000/2048));
                                    circlin_corr = zeros(height(angles_subwin),2);
                                    for t = 1:height(angles_subwin)
                                        [cc, pv] = circ_corrcl(1:decrease_count,angles_subwin(t,idx-decrease_count:idx-1));
                                        circlin_corr(t,1) = cc; circlin_corr(t,2) = pv;
                                    end
                                    %sum_sq_res = 0;
                                    %for ct = ct_decr(1):ct_decr(2)
                                    %    sum_sq_res = sum_sq_res + (peak_times(ct) - (coefs_decr(1)*ct + coefs_decr(2))).^2;
                                    % end
                                    traveling_ERP_table(end + 1, :) = {alignment, "peak", subwin_start, subwin_end, [idx-decrease_count idx-1], -1, peak_times, decr_speed, circlin_corr, ERP_table{row,1}};
                                    % {alignment, "peak", subwin_start, subwin_end, ct_decr, -1, peak_times, decr_speed, sum_sq_res_decr, circlin_corr_decr};
                                    if (cnum >= idx-decrease_count) & (cnum <= idx-1) % TEMP added 10/10/25
                                        travel_ERP_count = travel_ERP_count + 1;
                                        per_contact_wave_counts(cnum) = per_contact_wave_counts(cnum) + 1;
                                        trav_pks_bool(row) = 1;
                                    end
                                end
                                % if decrease_count > max_decr_count
                                %     max_decr_count = decrease_count; 
                                %     max_decr_end = idx; % max_decr_start = max_decr_end - max_decr_count
                                % end
                                decrease_count = 1; % reset
                            
                                while idx <= numel(peak_times) && peak_times(idx) == prev_latency
                                    prev_latency = peak_times(idx);
                                    idx = idx + 1;
                                end
                            
                            end
                            
                            % OG way of computing speed
                            % if max_incr_count >= 4
                            %     incr_speed = (max_incr_count - 1) / ((peak_times(max_incr_end - 1) - peak_times(max_incr_end - max_incr_count)) * (1000/2048));
                            % end
                            % if max_decr_count >= 4
                            %     decr_speed = (max_decr_count - 1) / ((peak_times(max_decr_end - 1) - peak_times(max_decr_end - max_decr_count)) * (1000/2048));
                            % end
                            
                            % if max_incr_count >= 4
                            %     ct_incr = [max_incr_end - max_incr_count, max_incr_end - 1];   
                            % end
                            % if max_decr_count >= 4
                            %     ct_decr = [max_decr_end - max_decr_count, max_decr_end - 1];
                            % end
    
                        end % end loop through rows in ERP_peaks_align
                    
                    end % end IF peaks non-empty
                %end % end loop through ref_phase

                if ~isempty(ERP_troughs_align)
                    ERP_troughs_align.IdxBefore = zeros(height(ERP_troughs_align),1);
                    ERP_troughs_align.IdxAfter = zeros(height(ERP_troughs_align),1);
                    for row = 1:height(ERP_troughs_align)
                        trough_idx = ERP_troughs_align{row,"Idx"};
                        name = ERP_troughs_align{row,"Contact"};
                        ct = find(strcmp(elec_name,name));
                        cnum = find(strcmp(cur_elec_contact_names,name));
                        per_contact_ERP_counts(cnum) = per_contact_ERP_counts(cnum) + 1;
                        
                        fname_data_avg = sprintf("%s/realavg_%s_%s_ch%03d.mat", data_out_dir, alignment, name, ct);
                        load(fname_data_avg, "avg_events")

                        fname_xgp_avg = sprintf("%s/centered_%s_xgp_%s_ch%03d.mat", data_out_dir, alignment, name, ct);
                        load(fname_xgp_avg, "xgp_avg")
    
                        % find all peaks (angle xgp method)
                        angle_posneg = xgp_avg;
                        angle_posneg(angle_posneg >= 0) = 1;
                        angle_posneg(angle_posneg < 0) = -1;
                        peaks_idx_align = find(diff(angle_posneg) > 0);
                        peaks_align = avg_events(peaks_idx_align);

                        % % find all peaks (deriv method)
                        % first_deriv_align = diff(avg_events);
                        % first_deriv_align(first_deriv_align < 0) = -1;
                        % first_deriv_align(first_deriv_align > 0) = 1;
                        % peaks_idx_align = find(diff(first_deriv_align) < 0);
                        % peaks_align = avg_events(peaks_idx_align);
    
                        % peaks before/after
                        idx_peak_before = peaks_idx_align(peaks_idx_align < trough_idx);
                        if ~isempty(idx_peak_before)
                            idx_peak_before = idx_peak_before(end);
                        else
                            idx_peak_before = 1;
                        end
                        ERP_troughs_align{row,"IdxBefore"} = idx_peak_before;
                        
                        idx_peak_after = peaks_idx_align(peaks_idx_align > trough_idx);
                        if ~isempty(idx_peak_after)
                            idx_peak_after = idx_peak_after(1);
                        else
                            idx_peak_after = numel(avg_events);
                        end
                        ERP_troughs_align{row,"IdxAfter"} = idx_peak_after;

                        % Test each row: take [IdxBefore, IdxAfter] as the subwindow
                        subwin_start = idx_peak_before;
                        subwin_end = idx_peak_after;
                        

                        xgp_avg_cur_elec = []; % each column = 1 contact on current electrode, each row = 1 timept
                        avg_cur_elec = [];
                        j = 1;
                        if ~contains(cur_letter, "'")
                            while ~contains(elec_name{j},cur_letter) || contains(elec_name{j},"'"), j = j+1; end
                        else
                            while ~contains(elec_name{j},cur_letter), j = j+1; end
                        end
                        cur_name = elec_name{j};
                        cur_elec_contact_names = []; % stores name of contacts on cur electrode
                        cur_elec_contact_ind = []; % stores indices of contacts on cur elec
                        
                        while strcmp(cur_name(1:2),cur_letter) || strcmp(cur_name(1),cur_letter)
                            cur_elec_contact_names = [cur_elec_contact_names; convertCharsToStrings(cur_name)];
                            cur_elec_contact_ind = [cur_elec_contact_ind; j];
                            
                            out_dir = sprintf('/media/Data/Human_Intracranial_MAD/analysis/Cross_Freq_Coupling/%s/%s/synchrony_data/generalized_phase',subject_ID, reference);
                            data_dir = sprintf('%s/data_filtered_trial_avg/%s',out_dir, alignment);
                            fname_xgp = sprintf("%s/centered_%s_xgp_%s_ch%03d.mat", data_dir, alignment, cur_name, j);
                            load(fname_xgp)
                            xgp_avg_cur_elec = [xgp_avg_cur_elec xgp_avg];
                    
                            fname_data_avg = sprintf("%s/realavg_%s_%s_ch%03d.mat", data_dir, alignment, cur_name, j);
                            load(fname_data_avg) % loads avg_events
                            %avg_cur_elec = [avg_cur_elec avg_events];
                        
                            j = j+1; 
                            if j > numel(elec_name), break; end
                            cur_name = elec_name{j};
                        end
                        angles_diff_sign_elec = angle(xgp_avg_cur_elec); % t rows, c columns (c = # of contacts)
                    
                       
                        angles_subwin = angles_diff_sign_elec(subwin_start:subwin_end,:);
                        %data_subwin = avg_cur_elec(subwin_starts(j):subwin_ends(j),:);
                
                        % FIND TIME of TROUGHs
                        ref_phase = pi;
                        peak_times = zeros(1,numel(cur_elec_contact_names)); % TROUGHS
                        for k = 1:numel(cur_elec_contact_names)
                            [peak1, peak_time] = min(abs(abs(angles_subwin(:,k)) - ref_phase));
                            % -pi is also acceptable, so angles_subwin(:,k) is wrapped in absolute value
                
                            peak_times(k) = peak_time;
                        end
                        
                        % FIND MONOTONE SEQUENCES
                        max_incr_count = 1; max_decr_count = 1;
                        max_incr_end = 1; max_decr_end = 1; %max_incr_st = 1; max_decr_st = 1;
                        increase_count = 1; decrease_count = 1;
                        prev_latency = peak_times(1);
                        idx = 2;
                        while idx <= numel(peak_times)
                            while idx <= numel(peak_times) && peak_times(idx) > prev_latency 
                                if ~(peak_times(idx) < subwin_end && prev_latency > 1)
                                    prev_latency = peak_times(idx); idx = idx + 1;
                                    break;
                                end
                                increase_count = increase_count + 1;
                                prev_latency = peak_times(idx);
                                idx = idx + 1;
                            end
                            if increase_count >= 4
                                %incr_wave_idx = [incr_wave_idx; [idx-increase_count idx-1]];
                                p = polyfit(1:increase_count,peak_times(idx-increase_count:idx-1),1);
                                incr_speed = 3.5 / (p(1) * (1000/2048));
                                circlin_corr = zeros(height(angles_subwin),2);
                                for t = 1:height(angles_subwin)
                                    [cc, pv] = circ_corrcl(1:increase_count,angles_subwin(t,idx-increase_count:idx-1));
                                    circlin_corr(t,1) = cc; circlin_corr(t,2) = pv;
                                end
                                traveling_ERP_table(end + 1, :) = {alignment, "trough", subwin_start, subwin_end, [idx-increase_count idx-1], 1, peak_times, incr_speed, circlin_corr, ERP_table{row,1}};
                                if (cnum >= idx-increase_count) & (cnum <= idx-1) % TEMP added 10/10/25
                                    travel_ERP_count = travel_ERP_count + 1;
                                    per_contact_wave_counts(cnum) = per_contact_wave_counts(cnum) + 1;
                                    trav_tghs_bool(row) = 1;
                                end
                            end
                            % if increase_count > max_incr_count
                            %     max_incr_count = increase_count; 
                            %     max_incr_end = idx;
                            % end
                            increase_count = 1; % reset
                           
                            while idx <= numel(peak_times) && peak_times(idx) < prev_latency 
                                if ~(peak_times(idx) > 1 && prev_latency < subwin_end)
                                    prev_latency = peak_times(idx); idx = idx + 1;
                                    break;
                                end
                                decrease_count = decrease_count + 1;
                                prev_latency = peak_times(idx);
                                idx = idx + 1;
                            end
                            if decrease_count >= 4
                                %decr_wave_idx = [decr_wave_idx; [idx-decrease_count idx-1]];
                                p = polyfit(1:decrease_count,peak_times(idx-decrease_count:idx-1),1);
                                decr_speed = 3.5 / (p(1) * (1000/2048));
                                circlin_corr = zeros(height(angles_subwin),2);
                                for t = 1:height(angles_subwin)
                                    [cc, pv] = circ_corrcl(1:decrease_count,angles_subwin(t,idx-decrease_count:idx-1));
                                    circlin_corr(t,1) = cc; circlin_corr(t,2) = pv;
                                end
                                traveling_ERP_table(end + 1, :) = {alignment, "trough", subwin_start, subwin_end, [idx-decrease_count idx-1], -1, peak_times, decr_speed, circlin_corr, ERP_table{row,1}};
                                if (cnum >= idx-decrease_count) & (cnum <= idx-1) % TEMP added 10/10/25
                                    travel_ERP_count = travel_ERP_count + 1;
                                    per_contact_wave_counts(cnum) = per_contact_wave_counts(cnum) + 1;
                                    trav_tghs_bool(row) = 1;
                                end
                            end
                            % if decrease_count > max_decr_count
                            %     max_decr_count = decrease_count; 
                            %     max_decr_end = idx; % max_decr_start = max_decr_end - max_decr_count
                            % end
                            decrease_count = 1; % reset
                        
                            while idx <= numel(peak_times) && peak_times(idx) == prev_latency
                                prev_latency = peak_times(idx);
                                idx = idx + 1;
                            end
                        
                        end
    
                    end % end loop through rows of trough ERP table
                    
                end % end IF trough table non-empty

                % % add trav_tghs_bool and trav_pks_bool onto ERP_peaks_align and ERP_troughs_align
                % ERP_peaks_align.TravWaveBool = trav_pks_bool;
                % ERP_troughs_align.TravWaveBool = trav_tghs_bool;

                % Compare (avg) amplitude for rows where trav_pks_bool = 1
                % with that for rows where trav_pks_bool = 0
                    % Commented out 10/30/25
                mean_amp_wave_pk = median(ERP_peaks_align{trav_pks_bool == 1,"Amplitude"});
                mean_amp_no_wave_pk = median(ERP_peaks_align{trav_pks_bool == 0,"Amplitude"});
                mean_amp_wave_tgh = median(ERP_troughs_align{trav_tghs_bool == 1,"Amplitude"});
                mean_amp_no_wave_tgh = mean(ERP_troughs_align{trav_tghs_bool == 0,"Amplitude"});
                fprintf("Avg amp of ERP peaks in waves: %f \n", mean_amp_wave_pk)
                fprintf("Avg amp of ERP peaks NOT in waves: %f \n", mean_amp_no_wave_pk)
                fprintf("Avg amp of ERP troughs in waves: %f \n", mean_amp_wave_tgh)
                fprintf("Avg amp of ERP troughs NOT in waves: %f \n", mean_amp_no_wave_tgh)

                total_ERPs = (height(ERP_troughs_align) + height(ERP_peaks_align));
                % %TEMP commented out 10/23/25
                %fprintf("%d / %d\n",travel_ERP_count, total_ERPs)

                % % TEMP commented out 10/30/25
                % for contact = 1:numel(per_contact_wave_counts)
                %     fprintf("%s: %d / %d \n", cur_elec_contact_names(contact), per_contact_wave_counts(contact), per_contact_ERP_counts(contact))
                % end
                % fprintf("*******************************************\n")

            end % end loop through alignments

            out_dir = sprintf('/media/Data/Human_Intracranial_MAD/analysis/Cross_Freq_Coupling/%s/%s/synchrony_data/generalized_phase',subject_ID, reference);
            table_out_dir = sprintf('%s/phase_progression_detection/trialavg',out_dir);
            if ~exist(table_out_dir, "dir"), mkdir(table_out_dir); end
            %traveling_ERP_table = cell2table(traveling_ERP_table,"VariableNames",["Alignment","Peak/Trough","Subwin St","Subwin End","Contacts","Direction","Latencies","Speed","SumSqRes","Corrs"]);
            traveling_ERP_table = cell2table(traveling_ERP_table,"VariableNames",["Alignment","Peak/Trough","Subwin St","Subwin End","Contacts","Direction","Latencies","Speed (m/s)","Corrs", "Ref Ct"]);
            fname_table = sprintf('%s/%s_trialavg_travERP_table.mat', table_out_dir, cur_letter); 
            save(fname_table, "traveling_ERP_table")
            %disp("traveling ERP table saved")
        end % end loop through elecs

    %end % end loop through alignments
    
end % end loop through subjects