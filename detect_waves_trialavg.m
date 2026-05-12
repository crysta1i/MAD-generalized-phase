%% TODO: create function that works for all references (Ground, neighbor_average)

% Input:
% - subwin_st, subwin_end
% - cur_letter ( = 'N')
% - phase (peak = 0, trough = pi)

cur_letter = 'N'; subwin_st = 200; subwin_end = 400; ref_phase = 0;
alignment = 'full_single_opt_info'; %fourth_unique_attribute
data_dir = sprintf('/media/Data/Human_Intracranial_MAD/analysis/Cross_Freq_Coupling/EMU038/neighbor_average/synchrony_data/generalized_phase/data_filtered_trial_avg/%s',alignment);
data_base_dir = '/media/Data/Human_Intracranial_MAD/1_formatted';

% Prerequisites
% - saved data_avg and xgp_avg for all contcats in cur_letter in trial_avg_traveling

if strcmp(reference,'Ground')
    load(sprintf('%s/%s/%s_MAD_SES%d_Setup.mat',data_base_dir,subject_ID,subject_ID,1));
    channel_ind = (1:numel(elec_name))';
else
    load(sprintf('%s/%s/%s_MAD_SES%d_Setup_%s.mat',data_base_dir,subject_ID,subject_ID,1,reference));
    elec_name = channel_name; elec_area = channel_area; elec_region = channel_region; elec_location = channel_location;
end

xgp_avg_cur_elec = []; % each column = 1 contact on current electrode, each row = 1 timept

addpath('/media/Data/Human_Intracranial_MAD/_toolbox');
[cur_elec_contact_ind, cur_elec_contact_names] = get_single_probe_contacts(reference, subject_ID, cur_letter);

for ct = 1:numel(cur_elec_contact_ind)
    contact = cur_elec_contact_ind(ct);
    cur_name = cur_elec_contact_names(ct);
    fname_data_avg = sprintf("%s/centered_%s_xgp_%s_ch%03d.mat", data_dir, alignment, cur_name, contact);
    load(fname_data_avg)
    xgp_avg_cur_elec = [xgp_avg_cur_elec xgp_avg];
end

angles_diff_sign_elec = angle(xgp_avg_cur_elec);
angles_subwin = angles_diff_sign_elec(subwin_st:subwin_end,:);

% FIND TIME of PEAKS
peak_times = zeros(1,numel(cur_elec_contact_names)); % PEAKS
for k = 1:numel(cur_elec_contact_names)
    [peak1, peak_time] = min(abs(abs(angles_subwin(:,k)) - ref_phase));
    peak_times(k) = peak_time;
end


%% 1. Peak detection on trial-averaged data

subject_IDs = {'EMU024','EMU025','EMU030','EMU038','EMU039','EMU041'};
reference = 'Ground'; sesnum = 1; Fs = 2048;
alignment = 'inspection'; % inspection % outcome % trialstart
    % amount % probability % novel_attribute_inspection % repeated_inspection
win = [0 1024]; %[0 1024] [1024 0]
datapath = '/media/Data/Human_Intracranial_MAD/';
data_base_dir = sprintf('%s1_formatted',datapath);
addpath('/media/Data/Human_Intracranial_MAD/analysis/PAC_code/CF_Coupling/generalized phase/generalized-phase/analysis');
addpath('/media/Data/Human_Intracranial_MAD/analysis/PAC_code/CF_Coupling/generalized phase/circstat-matlab')

for subject_num = 4 %:numel(subject_IDs) % loop through subjects
    subject_ID = subject_IDs{subject_num};

    % all setup files
    final_out_dir = sprintf('/media/Data/Human_Intracranial_MAD/analysis/Cross_Freq_Coupling/%s/%s/synchrony_data/generalized_phase',subject_ID, reference);
    table_out_dir = sprintf('%s/phase_progression_detection/trialavg',final_out_dir);
    if ~exist(table_out_dir, 'dir'), mkdir(table_out_dir); end

    out_dir = sprintf('/media/Data/Human_Intracranial_MAD/analysis/Cross_Freq_Coupling/%s/%s/synchrony_data',subject_ID, reference);  
    data_dir = sprintf('%s/generalized_phase/data_filtered_trial_avg/%s',out_dir, alignment);

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


    phase_prog_table = cell(0,8);

    for elec_num = 1:numel(elec_letters) % loop through electrodes
        
        cur_letter = elec_letters(elec_num);

        xgp_avg_cur_elec = []; % each column = 1 contact on current electrode, each row = 1 timept
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
            
            fname_data_avg = sprintf("%s/centered_%s_xgp_%s_ch%03d.mat", data_dir, alignment, cur_name, j);
            load(fname_data_avg)
            xgp_avg_cur_elec = [xgp_avg_cur_elec xgp_avg];

            j = j+1; 
            if j > numel(elec_name), break; end
            cur_name = elec_name{j};
        end
        angles_diff_sign_elec = angle(xgp_avg_cur_elec); % t rows, c columns (c = # of contacts)
        
        % EMU030 feedback: positive potential @350-500 samples
        % EMU038 inspection: F (+,225-375), L (+,400-500), M (+,200-300;
        % -,350-500), N (+, 200-325), P' (400-525), Y (400-500)

        %for subwin_start = 350 
            %subwindow = [subwin_start subwin_start + 150];
            %angles_cur_elec = angle(xgp_avg_cur_elec);
            %angles_cur_elec = angles_cur_elec(subwindow(1):subwindow(2),:);
            for phase_type = [0 1] % check both peak 0 and trough pi at each subwindow
                ref_phase = pi*phase_type;

                % TEMP %
                if phase_type == 1 % troughs come after peaks in EMU030 outcome
                    peak_trough = "trough";
                    subwin_start = 450;
                    subwindow = [450 650]; % for all elec EMU030 outcome
                    %subwindow = [450 600]; % for elec F,L,H',M',Z' EMU038 outcome

                    if strcmp(subject_ID, 'EMU038')
                        % FOR OUTCOME:
                        % if strcmp('W',cur_letter)
                        %     subwin_start = 400;
                        %     subwindow = [400 550];
                        % elseif strcmp('X',cur_letter)
                        %     subwin_start = 350;
                        %     subwindow = [350 550];
                        % elseif strcmp('Z',cur_letter)
                        %     subwin_start = 450;
                        %     subwindow = [450 600];
                        % end

                        % FOR INSPECTION:
                        if strcmp('M',cur_letter)
                            subwin_start = 350;
                            subwindow = [350 500];
                        else
                            continue;
                        end
                    end
                else
                    peak_trough = "peak";
                    subwin_start = 350;
                    subwindow = [350 500]; % for all elec EMU030 outcome
                    %subwindow = [350 500]; % for elec F,L,H',M',Z' EMU038 outcome

                    if strcmp(subject_ID, 'EMU038')
                        % FOR OUTCOME:
                        % if strcmp('W',cur_letter)
                        %     subwin_start = 300;
                        %     subwindow = [300 450];
                        % elseif strcmp('X',cur_letter)
                        %     subwin_start = 420;
                        %     subwindow = [420 600];
                        % elseif strcmp('Z',cur_letter)
                        %     continue; % only a negative potential (trough) visible on this elec
                        % end

                        % FOR INSPECTION:
                        if strcmp('L',cur_letter)
                            subwin_start = 400;
                            subwindow = [400 500];
                        elseif strcmp('M',cur_letter)
                            subwin_start = 200;
                            subwindow = [200 300];
                        elseif strcmp('F',cur_letter)
                            subwin_start = 225;
                            subwindow = [225 375];
                        elseif strcmp('N',cur_letter)
                            subwin_start = 200;
                            subwindow = [200 325];
                        elseif strcmp("P'",cur_letter)
                            subwin_start = 400;
                            subwindow = [400 525];
                        elseif strcmp('Y',cur_letter)
                            subwin_start = 400;
                            subwindow = [400 500];
                        else
                            break;
                        end
                    end
                end
                angles_subwin = angles_diff_sign_elec(subwindow(1):subwindow(2),:);

                peak_times = zeros(1,width(angles_subwin));
                for e = 1:width(angles_subwin)
                    [peak1, peak_time] = min(abs(angles_subwin(:,e) - ref_phase));
                    peak_times(e) = peak_time;
                end
                
                
                % 1. Detect monotone sequences
                    % - Take diff 
                    % - Mark all entries >0 with 1 and all entries <0 with -1 
                    % - Determine if [-1 -1 -1 -1] and/or [1 1 1 1] are sublists of the resulting list
                diff_latencies = diff(peak_times);
                diff_latencies(diff_latencies > 0) = 1;
                
                    % stricter criteria: constrain gaps between peak times
                for d = 1:numel(diff_latencies)
                    if d == 1
                        diff_latencies(d) = 1;
                    end
                    if diff_latencies > 0
                        if abs(diff_latencies(d) - diff_latenceis(d-1)) < 15
                            diff_latencies(d) = 1;
                            diff_latencies(d-1) = 1;
                        else
                            diff_latencies(d) = 0;
                        end
                    end
                end

                diff_latencies(diff_latencies < 0) = -5;
                diff_latencies = [0 diff_latencies];
                consec_decr = -1*ones(1,3); consec_incr = ones(1,3);
                is_traveling_decr = contains(num2str(diff_latencies),num2str(consec_decr));
                is_traveling_incr = contains(num2str(diff_latencies),num2str(consec_incr));

                % Three 1s or -1s in a row indicates four contacts
                % participating in a traveling wave

                % Ex:
                % peak_times = [9 20 25 33 30 34 39 42 32 25 23 23 31]
                % diff_latencies = [11 5 8 -3 4 5 3 -10 -7 -2 0 8]
                % diff_latencies = [1 1 1 -1 1 1 1 -1 -1 -1 0 1] 
                % find(diff_latencies == -1) = [4 8 9 10]
                % find(diff_latencies == 1) = [1 2 3 5 6 7 12]

                if is_traveling_incr && is_traveling_decr
                    direction = 0;
                elseif is_traveling_incr
                    ct_decr = find(diff_latencies == 1);
                    direction = 1;
                elseif is_traveling_decr
                    ct_incr = find(diff_latencies == -1);
                    direction = -1;
                end

                % 2. Test for linearity: Phase angle vs distance at each time point
                                            % Phase latency vs distance... not well defined
                
                subwin_duration = subwindow(2)-subwindow(1)+1;
                phase_dist_corrs = zeros(subwin_duration,2);
                for time = 1:subwin_duration
                    phase_dist_corr = corrcoef(1:numel(peak_times),angles_subwin(time,:));
                    phase_dist_corr = phase_dist_corr(2,1);
                    phase_dist_corrs(time,1) = phase_dist_corr;

                    % Permutation testing at each time point??
                    corr_distr = zeros(5000,1); % distribution of corrcoefs
                    for perm = 1:5000
                        contactnums = 1:numel(peak_times);
                        shuffled = contactnums(randperm(length(contactnums)));
                        cur_perm = corrcoef(shuffled, angles_subwin(time,:));
                        corr_distr(perm) = abs(cur_perm(2,1));
                    end
                    corr_distr = sort(corr_distr,'descend');
                    obs = find(corr_distr <= abs(phase_dist_corr), 1, 'first');
                    if isempty(obs)
                        pvalue = 1; 
                    else
                        pvalue = (obs - 1) / 5000;
                    end
                    phase_dist_corrs(time,2) = pvalue < 0.05;
                end
                corr_sig = phase_dist_corrs;
                num_sig_timepts = nnz(corr_sig(:,2));

                % 3. Save results
                

                % % detect monotone sequences in this list
                % max_incr_count = 1; max_decr_count = 1;
                % max_incr_end = 1; max_decr_end = 1;
                % increase_count = 1; decrease_count = 1;
                % prev_latency = peak_times(1);
                % prev_diff = abs(peak_times(1) - peak_times(2));
                % idx = 2;
                % while idx <= numel(peak_times)
                %     while idx <= numel(peak_times) && peak_times(idx) > prev_latency %&& abs(abs(peak_times(idx) - prev_latency) - prev_diff) < 15
                %         if ~(peak_times(idx) < 100 && prev_latency > 1)
                %             prev_latency = peak_times(idx); idx = idx + 1;
                %             break;
                %         end
                %         increase_count = increase_count + 1;
                %         prev_diff = abs(peak_times(idx) - prev_latency);
                %         prev_latency = peak_times(idx);
                %         idx = idx + 1;
                %     end
                %     if increase_count > max_incr_count
                %         max_incr_count = increase_count; 
                %         max_incr_end = idx;
                %     end
                %     increase_count = 1; % reset
                %     %prev_diff = abs(peak_times(idx) - peak_times(idx + 1));
                % 
                %     while idx <= numel(peak_times) && peak_times(idx) < prev_latency %&& abs(abs(prev_latency - peak_times(idx)) - prev_diff) < 15
                %         if ~(peak_times(idx) > 1 && prev_latency < 100)
                %             prev_latency = peak_times(idx); idx = idx + 1;
                %             break;
                %         end
                %         decrease_count = decrease_count + 1;
                %         %prev_diff = abs(peak_times(idx) - prev_latency);
                %         prev_latency = peak_times(idx);
                %         idx = idx + 1;
                %     end
                %     if decrease_count > max_decr_count
                %         max_decr_count = decrease_count; 
                %         max_decr_end = idx; % max_decr_start = max_decr_end - max_decr_count
                %     end
                %     decrease_count = 1; % reset
                %     % prev_diff = abs(peak_times(idx) - peak_times(idx + 1));
                % 
                %     while idx <= numel(peak_times) && peak_times(idx) == prev_latency
                %         prev_latency = peak_times(idx);
                %         idx = idx + 1;
                %     end
                % 
                % end


                % if max_incr_count >= 4
                %     incr_speed = (max_incr_count - 1) / ((peak_times(max_incr_end - 1) - peak_times(max_incr_end - max_incr_count)) * (1000/2048));
                % end
                % if max_decr_count >= 4
                %     decr_speed = (max_decr_count - 1) / ((peak_times(max_decr_end - 1) - peak_times(max_decr_end - max_decr_count)) * (1000/2048));
                % end

                % Ensure that peak timings follow roughly linear progression:

                % incr_cond = max_incr_count >= 4; %&& islinear(); %range(diff(peak_times((max_incr_end - max_incr_count):(max_incr_end - 1)))) <= 20;
                % decr_cond = max_decr_count >= 4; %&& islinear(); %range(diff(peak_times((max_decr_end - max_decr_count):(max_decr_end - 1)))) <= 20;
                
                if is_traveling_incr || is_traveling_decr 
                    % if incr_cond && decr_cond
                    %     direction = 0; % both increasing and decreasing latencies
                    %     ct_incr = [max_incr_end - max_incr_count, max_incr_end - 1];
                    %     ct_decr = [max_decr_end - max_decr_count, max_decr_end - 1];
                    %     if max_incr_count >= max_decr_count
                    %         speed = incr_speed;
                    %     else
                    %         speed = decr_speed;
                    %     end
                    % elseif incr_cond
                    %     ct_incr = [max_incr_end - max_incr_count, max_incr_end - 1];
                    %     ct_decr = [];
                    %     direction = 1; % increasing latencies
                    %     speed = incr_speed;
                    % else
                    %     ct_decr = [max_decr_end - max_decr_count, max_decr_end - 1];
                    %     ct_incr = [];
                    %     direction = -1; % decreasing latencies
                    %     speed = decr_speed;
                    % end
                    % % create corr_sig
                    % corr_sig = zeros(100,2); % zeros(33,2) TODO
                    % for kk = 1:100 % length of subwindow % 3:3:100
                    %     pl = angles_cur_elec(kk,:); % j is time point index
                    %     [ cc, pv ] = circ_corrcl(pl, 1:numel(cur_elec_contact_names)); % all in a row so exact coordinates don't matter
                    %     corr_sig(kk,:) = [cc, pv < 0.05];
                    % end
                    % % TODO: approximate linearity check before including in table
                    phase_prog_table(end + 1, :) = {cur_letter, peak_trough, ct_incr, ct_decr, direction, corr_sig, num_sig_timepts, peak_times};
                end
            end % end loop through phase type
        %end % end loop through subwindows
    
    end % end loop through electrodes
    phase_prog_table = cell2table(phase_prog_table,"VariableNames",["Electrode", "Phase Type", "Incr Contacts","Decr Contacts","Direction","Corr","Num Sig Timepts","Latencies"]);

    fname_table = sprintf('%s/%s_travwave_table.mat', table_out_dir, alignment); 
    save(fname_table, "phase_prog_table")
    disp("saved")
end % end loop through subjects

%% 2.1 Peak detection on trial-averaged data based on significant ERPs
% Subwindows in which we look for traveling waves depend on time of
% peak/troughs before and after the significant ERP

reference = 'Ground';
alignments = ["inspection","single_opt_first_inspection","full_single_opt_info"];
datapath = '/media/Data/Human_Intracranial_MAD/';
data_base_dir = sprintf('%s1_formatted',datapath);
addpath('/media/Data/Human_Intracranial_MAD/analysis/PAC_code/CF_Coupling/generalized phase/circstat-matlab')

subject_IDs = {'EMU001','EMU024','EMU025','EMU030','EMU037','EMU038','EMU039','EMU041','EMU047','EMU051'};

for subject_num = 6 %1:numel(subject_IDs)
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
        
        for e = 1:numel(elec_letters)

            cur_letter = elec_letters(e);

            load(sprintf('%s/fdr_ERP_sig_table_%s.mat', out_dir, cur_letter))

            ERP_peaks = ERP_table(ERP_table{:,"Peak/Trough"} == "peak",:);
            ERP_troughs = ERP_table(ERP_table{:,"Peak/Trough"} == "trough",:);
            % ERP_peaks = sortrows(ERP_peaks,4,"ascend"); % sort Idx column
            % ERP_troughs = sortrows(ERP_troughs,4,"ascend"); % sort Idx column

            % Create new table for saving traveling waves detected
            traveling_ERP_table = cell(0,9);
            

            for alignment = alignments

                travel_ERP_count = 0;
                %fprintf("--------Proportion of ERPs in traveling waves--------\n")
                %fprintf("%s %s: ", alignment, cur_letter)

                out_dir = sprintf('/media/Data/Human_Intracranial_MAD/analysis/Cross_Freq_Coupling/%s/%s/synchrony_data/generalized_phase',subject_ID, reference);  
                final_out_dir = sprintf('%s/waveform_plots/trial_averaged/%s',out_dir, alignment);
                data_out_dir = sprintf('%s/data_filtered_trial_avg/%s',out_dir, alignment);
            
                ERP_peaks_align = ERP_peaks(strcmp(ERP_peaks{:,"Alignment"},alignment),:);
                ERP_troughs_align = ERP_troughs(strcmp(ERP_troughs{:,"Alignment"},alignment),:);
                
                % Delete all rows in ERP_peaks_align and ERP_troughs_align
                % that are not significant
                ERP_peaks_align(ERP_peaks_align{:,"Sig"}==0,:) = [];
                ERP_troughs_align(ERP_troughs_align{:,"Sig"}==0,:) = [];

                %if isempty(ERP_peaks_align) && isempty(ERP_troughs_align), continue; end

                    if ~isempty(ERP_peaks_align)
                        ERP_peaks_align.IdxBefore = zeros(height(ERP_peaks_align),1);
                        ERP_peaks_align.IdxAfter = zeros(height(ERP_peaks_align),1);
                        for row = 1:height(ERP_peaks_align)
                            peak_idx = ERP_peaks_align{row,"Idx"};
                            name = ERP_peaks_align{row,"Contact"};
                            ct = find(strcmp(elec_name,name));
                            
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
                            ct_ind = find(cur_elec_contact_ind == ct);

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
                                    traveling_ERP_table(end + 1, :) = {alignment, "peak", subwin_start, subwin_end, [idx-increase_count idx-1], 1, peak_times, incr_speed, circlin_corr};
                                    if (ct_ind >= idx-increase_count) & (ct_ind <= idx-1) % TEMP added 10/10/25
                                        travel_ERP_count = travel_ERP_count + 1;
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
                                    traveling_ERP_table(end + 1, :) = {alignment, "peak", subwin_start, subwin_end, [idx-decrease_count idx-1], -1, peak_times, decr_speed, circlin_corr};
                                    % {alignment, "peak", subwin_start, subwin_end, ct_decr, -1, peak_times, decr_speed, sum_sq_res_decr, circlin_corr_decr};
                                    if (ct_ind >= idx-decrease_count) & (ct_ind <= idx-1) % TEMP added 10/10/25
                                        travel_ERP_count = travel_ERP_count + 1;
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
                        cur_name = ERP_troughs_align{row,"Contact"};
                        contact = find(strcmp(elec_name,cur_name));
                        
                        fname_data_avg = sprintf("%s/realavg_%s_%s_ch%03d.mat", data_out_dir, alignment, cur_name, contact);
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
                                traveling_ERP_table(end + 1, :) = {alignment, "trough", subwin_start, subwin_end, [idx-increase_count idx-1], 1, peak_times, incr_speed, circlin_corr};
                                if (ct_ind >= idx-increase_count) & (ct_ind <= idx-1) % TEMP added 10/10/25
                                    travel_ERP_count = travel_ERP_count + 1;
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
                                traveling_ERP_table(end + 1, :) = {alignment, "trough", subwin_start, subwin_end, [idx-decrease_count idx-1], -1, peak_times, decr_speed, circlin_corr};
                                if (ct_ind >= idx-decrease_count) & (ct_ind <= idx-1) % TEMP added 10/10/25
                                    travel_ERP_count = travel_ERP_count + 1;
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
                
                total_ERPs = (height(ERP_troughs_align) + height(ERP_peaks_align));
                %fprintf("%d / %d\n",travel_ERP_count, total_ERPs)

            end % end loop through alignments

            % TEMP: Commented out on 10/9/25
            out_dir = sprintf('/media/Data/Human_Intracranial_MAD/analysis/Cross_Freq_Coupling/%s/%s/synchrony_data/generalized_phase',subject_ID, reference);
            table_out_dir = sprintf('%s/phase_progression_detection/trialavg',out_dir);
            %traveling_ERP_table = cell2table(traveling_ERP_table,"VariableNames",["Alignment","Peak/Trough","Subwin St","Subwin End","Contacts","Direction","Latencies","Speed","SumSqRes","Corrs"]);
            traveling_ERP_table = cell2table(traveling_ERP_table,"VariableNames",["Alignment","Peak/Trough","Subwin St","Subwin End","Contacts","Direction","Latencies","Speed (m/s)","Corrs"]);
            fname_table = sprintf('%s/%s_trialavg_travERP_table.mat', table_out_dir, cur_letter); 
            save(fname_table, "traveling_ERP_table")
        end % end loop through elecs

    %end % end loop through alignments
    
end % end loop through subjects

%% 2.2 Peak detection on trial-averaged data based on significant ERPs
% look for traveling waves, EXCLUDING DUPLICATES

reference = 'Ground';
alignments = ["inspection","single_opt_first_inspection","full_single_opt_info"];
datapath = '/media/Data/Human_Intracranial_MAD/';
data_base_dir = sprintf('%s1_formatted',datapath);
addpath('/media/Data/Human_Intracranial_MAD/analysis/PAC_code/CF_Coupling/generalized phase/circstat-matlab')

subject_IDs = {'EMU001','EMU024','EMU025','EMU030','EMU037','EMU038','EMU039','EMU041','EMU047','EMU051'};

for subject_num = 6 %1:numel(subject_IDs)
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
        
        for e = 1:numel(elec_letters)

            cur_letter = elec_letters(e);

            load(sprintf('%s/fdr_ERP_sig_table_%s.mat', out_dir, cur_letter))

            ERP_peaks = ERP_table(ERP_table{:,"Peak/Trough"} == "peak",:);
            ERP_troughs = ERP_table(ERP_table{:,"Peak/Trough"} == "trough",:);
            % ERP_peaks = sortrows(ERP_peaks,4,"ascend"); % sort Idx column
            % ERP_troughs = sortrows(ERP_troughs,4,"ascend"); % sort Idx column

            % Create new table for saving traveling waves detected
            unique_trav_ERP_table = cell(0,9);
            

            for alignment = alignments

                travel_ERP_count = 0;
                %fprintf("--------Proportion of ERPs in traveling waves--------\n")
                %fprintf("%s %s: ", alignment, cur_letter)

                out_dir = sprintf('/media/Data/Human_Intracranial_MAD/analysis/Cross_Freq_Coupling/%s/%s/synchrony_data/generalized_phase',subject_ID, reference);  
                final_out_dir = sprintf('%s/waveform_plots/trial_averaged/%s',out_dir, alignment);
                data_out_dir = sprintf('%s/data_filtered_trial_avg/%s',out_dir, alignment);
            
                ERP_peaks_align = ERP_peaks(strcmp(ERP_peaks{:,"Alignment"},alignment),:);
                ERP_troughs_align = ERP_troughs(strcmp(ERP_troughs{:,"Alignment"},alignment),:);
                
                % Delete all rows in ERP_peaks_align and ERP_troughs_align
                % that are not significant
                ERP_peaks_align(ERP_peaks_align{:,"Sig"}==0,:) = [];
                ERP_troughs_align(ERP_troughs_align{:,"Sig"}==0,:) = [];

                %if isempty(ERP_peaks_align) && isempty(ERP_troughs_align), continue; end

                    if ~isempty(ERP_peaks_align)
                        ERP_peaks_align.IdxBefore = zeros(height(ERP_peaks_align),1);
                        ERP_peaks_align.IdxAfter = zeros(height(ERP_peaks_align),1);
                        for row = 1:height(ERP_peaks_align)
                            peak_idx = ERP_peaks_align{row,"Idx"};
                            name = ERP_peaks_align{row,"Contact"};
                            ct = find(strcmp(elec_name,name));
                            
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
                            ct_ind = find(cur_elec_contact_ind == ct);

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

                                    recounted = 0;
                                    for wave_ct = idx-increase_count:idx-1
                                        wave_ct_name = cur_elec_contact_names(wave_ct);
                                        wave_ct_idx = subwin_start + peak_times(wave_ct);
                                        filt_idx = strcmp(ERP_peaks_align.Contact,wave_ct_name) & ERP_peaks_align.Idx == wave_ct_idx;
                                        filt_idx(filt_idx >= row) = []; % only appearances in previous row indicate double-counting
                                        if ~all(filt_idx==0)
                                            recounted = 1;
                                            break;
                                        end
                                    end
                                    
                                    if recounted == 0
                                        unique_trav_ERP_table(end + 1, :) = {alignment, "peak", subwin_start, subwin_end, [idx-increase_count idx-1], 1, peak_times, incr_speed, circlin_corr};
                                        if (ct_ind >= idx-increase_count) & (ct_ind <= idx-1) % TEMP added 10/10/25
                                            travel_ERP_count = travel_ERP_count + 1;
                                        end
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

                                    
                                    recounted = 0;
                                    for wave_ct = idx-decrease_count:idx-1
                                        wave_ct_name = cur_elec_contact_names(wave_ct);
                                        wave_ct_idx = subwin_start + peak_times(wave_ct);
                                        filt_idx = strcmp(ERP_peaks_align.Contact,wave_ct_name) & ERP_peaks_align.Idx == wave_ct_idx;
                                        filt_idx(filt_idx >= row) = []; 
                                        if ~all(filt_idx==0)
                                            recounted = 1;
                                            break;
                                        end
                                    end
                                    
                                    if recounted == 0
                                        unique_trav_ERP_table(end + 1, :) = {alignment, "peak", subwin_start, subwin_end, [idx-decrease_count idx-1], -1, peak_times, decr_speed, circlin_corr};
                                        if (ct_ind >= idx-decrease_count) & (ct_ind <= idx-1) % TEMP added 10/10/25
                                            travel_ERP_count = travel_ERP_count + 1;
                                        end
                                    end

                                end

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
                        cur_name = ERP_troughs_align{row,"Contact"};
                        contact = find(strcmp(elec_name,cur_name));
                        
                        fname_data_avg = sprintf("%s/realavg_%s_%s_ch%03d.mat", data_out_dir, alignment, cur_name, contact);
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

                                recounted = 0;
                                for wave_ct = idx-increase_count:idx-1
                                    wave_ct_name = cur_elec_contact_names(wave_ct);
                                    wave_ct_idx = subwin_start + peak_times(wave_ct);
                                    filt_idx = strcmp(ERP_troughs_align.Contact,wave_ct_name) & ERP_troughs_align.Idx == wave_ct_idx;
                                    filt_idx(filt_idx >= row) = []; 
                                    if ~all(filt_idx==0)
                                        recounted = 1;
                                        break;
                                    end
                                end
                                
                                if recounted == 0
                                    unique_trav_ERP_table(end + 1, :) = {alignment, "trough", subwin_start, subwin_end, [idx-increase_count idx-1], 1, peak_times, incr_speed, circlin_corr};
                                    if (ct_ind >= idx-increase_count) & (ct_ind <= idx-1) % TEMP added 10/10/25
                                        travel_ERP_count = travel_ERP_count + 1;
                                    end
                                end
                                
                            end
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
                                
                                recounted = 0;
                                for wave_ct = idx-decrease_count:idx-1
                                    wave_ct_name = cur_elec_contact_names(wave_ct);
                                    wave_ct_idx = subwin_start + peak_times(wave_ct);
                                    filt_idx = strcmp(ERP_troughs_align.Contact,wave_ct_name) & ERP_troughs_align.Idx == wave_ct_idx;
                                    filt_idx(filt_idx >= row) = []; 
                                    if ~all(filt_idx==0)
                                        recounted = 1;
                                        break;
                                    end
                                end
                                
                                if recounted == 0
                                    unique_trav_ERP_table(end + 1, :) = {alignment, "trough", subwin_start, subwin_end, [idx-decrease_count idx-1], -1, peak_times, decr_speed, circlin_corr};
                                    if (ct_ind >= idx-decrease_count) & (ct_ind <= idx-1) % TEMP added 10/10/25
                                        travel_ERP_count = travel_ERP_count + 1;
                                    end
                                end
                            end
                            decrease_count = 1; % reset
                        
                            while idx <= numel(peak_times) && peak_times(idx) == prev_latency
                                prev_latency = peak_times(idx);
                                idx = idx + 1;
                            end
                        
                        end
    
                    end % end loop through rows of trough ERP table

                    
                end % end IF trough table non-empty
                
                total_ERPs = (height(ERP_troughs_align) + height(ERP_peaks_align));
                %fprintf("%d / %d\n",travel_ERP_count, total_ERPs)

            end % end loop through alignments

            out_dir = sprintf('/media/Data/Human_Intracranial_MAD/analysis/Cross_Freq_Coupling/%s/%s/synchrony_data/generalized_phase',subject_ID, reference);
            table_out_dir = sprintf('%s/phase_progression_detection/trialavg',out_dir);
            %traveling_ERP_table = cell2table(traveling_ERP_table,"VariableNames",["Alignment","Peak/Trough","Subwin St","Subwin End","Contacts","Direction","Latencies","Speed","SumSqRes","Corrs"]);
            unique_trav_ERP_table = cell2table(unique_trav_ERP_table,"VariableNames",["Alignment","Peak/Trough","Subwin St","Subwin End","Contacts","Direction","Latencies","Speed (m/s)","Corrs"]);
            fname_table = sprintf('%s/%s_trialavg_unqiue_travERPtable.mat', table_out_dir, cur_letter); 
            save(fname_table, "unique_trav_ERP_table")
        end % end loop through elecs

    %end % end loop through alignments
    
end % end loop through subjects

%% 3. Based on times of trial-avg traveling waves, look for traveling waves in single trial data
reference = 'Ground';
alignments = ["inspection","single_opt_first_inspection","full_single_opt_info"];
datapath = '/media/Data/Human_Intracranial_MAD/';
data_base_dir = sprintf('%s1_formatted',datapath);
addpath('/media/Data/Human_Intracranial_MAD/analysis/PAC_code/CF_Coupling/generalized phase/circstat-matlab')
addpath('/media/Data/Human_Intracranial_MAD/_toolbox/spectral-analysis')

subject_IDs = {'EMU001','EMU024','EMU025','EMU030','EMU037','EMU038','EMU039','EMU041','EMU047','EMU051'};

for subject_num = 6 %1:numel(subject_IDs)
    subject_ID = subject_IDs{subject_num}; 
    
    if strcmp(subject_ID,'EMU001')
        Fs = 1000;
    else
        Fs = 2048;
    end
    
    out_dir = sprintf('/media/Data/Human_Intracranial_MAD/analysis/Cross_Freq_Coupling/%s/%s/synchrony_data/generalized_phase',subject_ID, reference);    
    table_out_dir = sprintf('%s/phase_progression_detection/trialavg',out_dir);
            
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
    
    %hist_single_trial_ct = zeros(1,10);
    parsave = @(fname, x) save(fname, 'x');
   for e = 1:numel(elec_letters)
        cur_letter = elec_letters(e);

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


        % fname_table = sprintf('%s/%s_trialavg_travERP_table.mat', table_out_dir, cur_letter); 
        % load(fname_table, "traveling_ERP_table")

        fname_table = sprintf('%s/%s_trialavg_unqiue_travERPtable.mat', table_out_dir, cur_letter); 
        traveling_ERP_table = load(fname_table, "unique_trav_ERP_table");

        %fprintf('----------Elec %s # single trials with matching waves----------\n',cur_letter)
        
        n_trav_waves = min(height(traveling_ERP_table),15); % 15 traveling waves

        % % get rows with top 15 unique amplitudes
        % traveling_ERP_table = sortrows(traveling_ERP_table, "MaxAmp","descend"); 
        % high_amp_rows = [1]; cur_row = 2;
        % while (numel(high_amp_rows) < n_trav_waves) && (cur_row < height(traveling_ERP_table))
        %     cur_amp = traveling_ERP_table{cur_row, "MaxAmp"};
        %     prev_amp = traveling_ERP_table{cur_row - 1, "MaxAmp"};
        %     if cur_amp ~= prev_amp
        %         high_amp_rows = [high_amp_rows; cur_row];
        %     end
        %     cur_row = cur_row + 1;
        % end

        matching_single_trials = zeros(height(traveling_ERP_table),1);
        %rand_rows = randsample(height(traveling_ERP_table),n_trav_waves); 
        
        for r = 1:height(traveling_ERP_table) %numel(rand_rows) %numel(high_amp_rows)
            row = r;
            %row = rand_rows(r);
            %row = high_amp_rows(r);
            alignment = traveling_ERP_table{row,1};
            win = [0 Fs/2];
            
            % START RUNNING HERE 1/24/26
            single_trial_count = 0; % count number of trials with matching traveling waves
            trials_per_ses = 108/num_sessions; % 108 is divisible by every possible number of sessions (1, 2, 3, 4, 9)
            
            for sesnum = 1:num_sessions
                load(sprintf('%s/%s/%s_MAD_SES%d_Setup.mat',data_base_dir,subject_ID,subject_ID,sesnum));
                
                [align_times,trial_numbers] = get_align_times(filters, trial_times, trial_words, alignment);
                remove_ind = isnan(align_times);
                align_times(isnan(align_times)) = []; trial_numbers(remove_ind) = [];
                align_times = round(align_times*Fs);

                rng(123);
                rand_events = randsample(numel(align_times),trials_per_ses);
                for iter = 1:trials_per_ses  % loop through events (test 100 single trials)
                    rand_event = rand_events(iter);    
                    %rand_event = iter;

                    %subwin_st = 200; subwin_end = 400;
                    subwin_st = traveling_ERP_table{row,3};
                    subwin_end = traveling_ERP_table{row,4};
     
                    filt_signal = []; % each column = signal of one contact
                    peak_times = zeros(numel(cur_elec_contact_ind),1);
                    for cnum = 1:numel(cur_elec_contact_names)
                        cur_name = cur_elec_contact_names(cnum);
                        contact = cur_elec_contact_ind(cnum);
                        
                        % load real-valued signal (separate_channel_files -> filter)
                        data = load(sprintf('%s/%s/separate_channel_files/%s_MAD_SES%d_ch%03d.mat',data_base_dir,subject_ID,subject_ID,sesnum,contact)); % for saving phase for ALL contacts
                        [b,a] = butter( 4, [5 40] ./ (Fs/2) ); data = filtfilt( b, a, data.data );
                        data_event = data(align_times(rand_event)-win(1):align_times(rand_event)+win(2)-1);
                        data_event = ( data_event - mean(data_event) ) ./ std(data_event);
                        
                        filt_signal = [filt_signal data_event];

                        xgp = load(sprintf('%s/extracted_phase/xgp_%s_ses%02d_ch%03d.mat',out_dir,cur_name,sesnum,contact));
                        xgp_event = xgp(align_times(rand_event)-win(1):align_times(rand_event)+win(2)-1);
                        xgp_event_subwin = xgp_event(subwin_st:subwin_end);
                        angle_subwin = angle(xgp_event_subwin);

                        if strcmp(traveling_ERP_table{row,"Peak/Trough"},"peak")
                            [~,peak_time] = min(abs(angle_subwin));
                        else % trough
                            [~,peak_time] = min(abs(abs(angle_subwin) - pi));
                        end
                        %[~,peak_time] = min(abs(angle_subwin)); % 1.24.26 - only dealing with peaks, temp commented out block above

                        peak_times(cnum) = peak_time;
                        
                    end % end loop through contacts on one probe
                    filt_signal = filt_signal(subwin_st:subwin_end,:);

                    % % plot filt_signal and display it
                    % colors = [
                    %     0.76 0.86 0.98;
                    %     0.72 0.83 0.96;
                    %     0.68 0.79 0.94;
                    %     0.64 0.76 0.92;
                    %     0.60 0.72 0.90;
                    %     0.56 0.69 0.88;
                    %     0.52 0.65 0.86;
                    %     0.48 0.62 0.84;
                    %     0.44 0.58 0.82;
                    %     0.40 0.55 0.80;
                    %     0.36 0.51 0.78;
                    %     0.32 0.48 0.76;
                    %     0.28 0.44 0.74;
                    %     0.24 0.41 0.72;
                    %     0.20 0.37 0.70;
                    %     0.16 0.34 0.68;
                    %     0.12 0.30 0.66;
                    %     0.06 0.26 0.63];
                    % fg1 = figure; 
                    % set( fg1, 'position', [ 88  100  1400  600 ] ) %[ 88  1593  1250  420 ]
                    % for cnum = 1:numel(cur_elec_contact_ind)
                    %     plot(filt_signal(:,cnum),'linewidth',2,'color',colors(cnum,:)) % -1024:-1
                    %     hold on
                    %     xlim([subwin_st,subwin_end]); 
                    %     %ylim([-0.5, 0.5]);
                    %     xlabel('samples from alignment'); 
                    %     ylabel('zscored amplitude');
                    %     yl1 = ylim;
                    %     title(sprintf("%s Elec %s (%d-%d)",alignment,cur_letter, subwin_st, subwin_end))
                    % end
                    % 
                    % close;

                    % check for traveling waves in peak_times
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
                            % check whether traveling wave is in same direction,
                            % involving same contacts (any overlap in contacts is acceptable)

                            % % 1.24.26 -- require that all (4+) contacts
                            % % participating in traveling wave in trial-avg
                            % % are invovled in the single-trial wave
                            %     % HARD-CODE: (TEMP - delete this after 1.24.26)
                            % trial_avg_cts = [2,5];
                            % trial_avg_direction = -1;

                            ct_incr = [idx-increase_count idx-1];
                            trial_avg_direction = traveling_ERP_table{row,"Direction"};
                            trial_avg_cts = traveling_ERP_table{row,"Contacts"};
                            overlap = (ct_incr(1) >= trial_avg_cts(1) && ct_incr(1) <= trial_avg_cts(2)) || (ct_incr(2) >= trial_avg_cts(1) && ct_incr(2) <= trial_avg_cts(2));
                            
                            %overlap = (ct_incr(1) <= trial_avg_cts(1)) && (ct_incr(2) >= trial_avg_cts(2));
                            if (trial_avg_direction == 1) & (overlap)
                                single_trial_count = single_trial_count + 1;
                            end
                        end
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
                            % trial_avg_direction = -1;
                            % trial_avg_cts = [2,5];

                            ct_decr = [idx-decrease_count idx-1];
                            trial_avg_direction = traveling_ERP_table{row,"Direction"};
                            trial_avg_cts = traveling_ERP_table{row,"Contacts"};
                            overlap = (ct_decr(1) >= trial_avg_cts(1) && ct_decr(1) <= trial_avg_cts(2)) || (ct_decr(2) >= trial_avg_cts(1) && ct_decr(2) <= trial_avg_cts(2));
                            %overlap = (trial_avg_cts(1) >= ct_decr(1)) && (trial_avg_cts(2) <= ct_decr(2));
                            if (trial_avg_direction == -1) & (overlap)
                                single_trial_count = single_trial_count + 1;
                            end
                        end
                        decrease_count = 1; % reset
                    
                        while idx <= numel(peak_times) && peak_times(idx) == prev_latency
                            prev_latency = peak_times(idx);
                            idx = idx + 1;
                        end
                    
                    end

                end % end loop through events

            end % end loop through sessions

            %fprintf("Elec N: %d single trials with waves that match main one in trial-avg", single_trial_count)
            % END RUN HERE 1/24/26
            
            subwin_st = traveling_ERP_table{row,3};
            subwin_end = traveling_ERP_table{row,4};

            % fprintf("Wave %d (%s %d-%d): %d\n", row, alignment, subwin_st, subwin_end, single_trial_count)
            % if single_trial_count > 90
            %     hist_single_trial_ct(10) = hist_single_trial_ct(10) + 1;
            % elseif single_trial_count > 80
            %     hist_single_trial_ct(9) = hist_single_trial_ct(9) + 1;
            % elseif single_trial_count > 70
            %     hist_single_trial_ct(8) = hist_single_trial_ct(8) + 1;
            % elseif single_trial_count > 60
            %     hist_single_trial_ct(7) = hist_single_trial_ct(7) + 1;
            % elseif single_trial_count > 50
            %     hist_single_trial_ct(6) = hist_single_trial_ct(6) + 1;
            % elseif single_trial_count > 40
            %     hist_single_trial_ct(5) = hist_single_trial_ct(5) + 1;
            % elseif single_trial_count > 30
            %     hist_single_trial_ct(4) = hist_single_trial_ct(4) + 1;
            % elseif single_trial_count > 20
            %     hist_single_trial_ct(3) = hist_single_trial_ct(3) + 1;
            % elseif single_trial_count > 10
            %     hist_single_trial_ct(2) = hist_single_trial_ct(2) + 1;
            % else
            %     hist_single_trial_ct(1) = hist_single_trial_ct(1) + 1;
            % end

            matching_single_trials(r) = single_trial_count;

        end % end loop through rows in trav wave table
    
        traveling_ERP_table.SingleTrialCt = matching_single_trials;
        
        %fname_table = sprintf('%s/%s_trialavg_unqiue_travERPtable.mat', table_out_dir, cur_letter); 
        %parsave(fname_table,unique_trav_ERP_table)

        fname_table = sprintf('%s/%s_trialavg_travERPtable.mat', table_out_dir, cur_letter); 
        save(fname_table,"traveling_ERP_table")

    end % end loop through elecs
    %fprintf("****# of waves with 10-90 percent match in single trials****")
    %disp(hist_single_trial_ct)

end % end loop though subjects

%% 1.25.26 - # single trials with matching waves (single proble - elec N)
% Hard-coded the participating contacts, require exact match -- all contacts overlap

% START RUNNING HERE 1/24/26
addpath('/media/Data/Human_Intracranial_MAD/_toolbox');
datapath = '/media/Data/Human_Intracranial_MAD/';
data_base_dir = sprintf('%s1_formatted',datapath); 
reference = 'neighbor_average'; %'neighbor_average'; % Ground
subject_ID = 'EMU038';
alignment = 'fourth_unique_attribute'; %'full_single_opt_info'; 

out_dir = sprintf('/media/Data/Human_Intracranial_MAD/analysis/Cross_Freq_Coupling/%s/%s/synchrony_data/generalized_phase',subject_ID, reference);    
data_avg_dir = sprintf('%s/data_filtered_trial_avg/%s',out_dir,alignment);

peak_wave = true; % peak_wave = false if looking for traveling troughs
cur_letter = 'N'; % make sure cur_elec_contact_ind/name are created
trial_avg_cts = [7,10]; % N3-N6
ref_ct = trial_avg_cts(1);
trial_avg_direction = -1;
subwin_st = 200; subwin_end = 400;
if strcmp(subject_ID,'EMU001')
    Fs = 1000;
else
    Fs = 2048;
end
win = [0 Fs/2];

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

[cur_elec_contact_ind, cur_elec_contact_names] = get_single_probe_contacts(reference, subject_ID, cur_letter);

% load trial-avg data to find time and ct of highest peak
trial_avg_all = []; 
for cnum = trial_avg_cts(1):trial_avg_cts(2)
    contact = cur_elec_contact_ind(cnum);
    cur_name = cur_elec_contact_names(cnum);
    fname_data_avg = sprintf("%s/realavg_%s_%s_ch%03d.mat", data_avg_dir, alignment, cur_name, contact);
    load(fname_data_avg, "avg_events") 
    trial_avg_all = [trial_avg_all avg_events];
end

if peak_wave
    % highest-amplitude peak
    [~, linearIndex] = max(trial_avg_all, [], "all", "linear");
    % [time_max_pk, ct_max_pk] = ...
    [time_max_pk, ct_max_pk] = ind2sub(size(trial_avg_all), linearIndex);
    ct_max_pk = ct_max_pk + trial_avg_cts(1) - 1; % convert column indices in to_plot_normal to indices in list of contacts on current probe
    phase = "Peak";
else % align troughs
    % find largest-amplitude trough
    [~, linearIndex] = min(trial_avg_all, [], "all", "linear");
    % time_max_pk
    [time_max_pk, ct_max_pk] = ind2sub(size(trial_avg_all), linearIndex);
    ct_max_pk = ct_max_pk + trial_avg_cts(1) - 1;
    phase = "Trough";
end

single_trial_count = 0; % count number of trials with matching traveling waves
trials_per_ses = 108/num_sessions; % 108 is divisible by every possible number of sessions (1, 2, 3, 4, 9)

for sesnum = 1:num_sessions
    if strcmp(reference,'Ground')
        load(sprintf('%s/%s/%s_MAD_SES%d_Setup.mat',data_base_dir,subject_ID,subject_ID,sesnum));
    else
        load(sprintf('%s/%s/%s_MAD_SES%d_Setup_%s.mat',data_base_dir,subject_ID,subject_ID,sesnum,reference));
    end
    
    [align_times,trial_numbers] = get_align_times(filters, trial_times, trial_words, alignment);
    remove_ind = isnan(align_times);
    align_times(isnan(align_times)) = []; trial_numbers(remove_ind) = [];
    align_times = round(align_times*Fs);

    rng(123);
    rand_events = randsample(numel(align_times),trials_per_ses);
    
    single_trials = [];
    for iter = 1:numel(align_times)   %trials_per_ses   (TEMP commented out)   % loop through events (test 100 single trials)
        %rand_event = rand_events(iter);    % TEMP commented out
        rand_event = iter;

        % subwin_st = traveling_ERP_table{row,3};
        % subwin_end = traveling_ERP_table{row,4};

        raw_signal = []; zscore_signal = [];
        filt_signal = []; % each column = signal of one contact
        peak_times = zeros(numel(cur_elec_contact_ind),1);
        for cnum = 1:numel(cur_elec_contact_ind)
            cur_name = cur_elec_contact_names(cnum);
            contact = cur_elec_contact_ind(cnum);
            
            % load real-valued signal (separate_channel_files -> filter)
            if strcmp(reference,'Ground')
                load(sprintf('%s/%s/separate_channel_files/%s_MAD_SES%d_ch%03d.mat',data_base_dir,subject_ID,subject_ID,sesnum,contact)); 
            elseif strcmp(reference,'neighbor_average')
                load(sprintf('%s/%s/separate_channel_files/%s/%s_MAD_SES%d_ch%03d.mat',data_base_dir,subject_ID,reference,subject_ID,sesnum,contact));
            end
            
            raw = data; raw_event = raw(align_times(rand_event)-win(1):align_times(rand_event)+win(2)-1);
            zscored_event = (raw_event - mean(raw_event)) ./ std(raw_event);
            [b,a] = butter( 4, [5 40] ./ (Fs/2) ); data = filtfilt( b, a, data );
            data_event = data(align_times(rand_event)-win(1):align_times(rand_event)+win(2)-1);
            data_event_padded = data(align_times(rand_event)-win(1)-Fs:align_times(rand_event)+win(2)-1+Fs);
                % ^ padding to account for edge artifacts caused by generalized_phase_vector function
            data_event = ( data_event - mean(data_event) ) ./ std(data_event);
            
            raw_signal = [raw_signal raw_event];
            zscore_signal = [zscore_signal zscored_event];
            filt_signal = [filt_signal data_event];
    
            % 2/21/26 -- no longer using angle(xgp) method to find peaks/troughs
            % 3/6/26 -- back to using angle(xgp) method for finding peaks, after correcting typo
            xgp_padded = generalized_phase_vector(data_event_padded,Fs,5);
            xgp_event = xgp_padded(1+Fs:1+Fs+win(2)-1);
            xgp_event_subwin = xgp_event(subwin_st:subwin_end);

            % xgp = load(sprintf('%s/extracted_phase/xgp_%s_ses%02d_ch%03d.mat',out_dir,cur_name,sesnum,contact));
            % xgp = xgp.xgp;
            % xgp_event = xgp(align_times(rand_event)-win(1):align_times(rand_event)+win(2)-1);
            % xgp_event_subwin = xgp_event(subwin_st:subwin_end);
            angle_subwin = angle(xgp_event_subwin);
            trial_angle_pn = angle_subwin;
            trial_angle_pn(trial_angle_pn >= 0) = 1; trial_angle_pn(trial_angle_pn < 0) = -1;

            %cur_ct_pk_times = [];
            if peak_wave
                %cur_ct_pk_times = find(sign(diff([data_event; data_event(end)]))<0&sign(diff([data_event(1); data_event]))>0);
                cur_ct_pk_times = find(diff(trial_angle_pn) > 0);
                if isempty(cur_ct_pk_times)
                    [~,cur_ct_pk_times] = min(abs(angle_subwin));
                end
            else
                %cur_ct_pk_times = find(sign(diff([data_event; data_event(end)]))>0&sign(diff([data_event(1); data_event]))<0);
                cur_ct_pk_times = find(diff(trial_angle_pn) < 0);
                if isempty(cur_ct_pk_times)
                    [~,cur_ct_pk_times] = min(abs(abs(angle_subwin) - pi));
                end
            end
    
            
            [~,idx_mindiff] = min(abs(cur_ct_pk_times - time_max_pk));
            time_closest_pk = cur_ct_pk_times(idx_mindiff);

            % if strcmp(traveling_ERP_table{row,"Peak/Trough"},"peak")
            %     [~,peak_time] = min(abs(angle_subwin));
            %     %[~, peak_time] = max(data_event);
            % else % trough
            %     [~,peak_time] = min(abs(abs(angle_subwin) - pi));
            %     %[~,peak_time] = min(data_event);
            % end
            
            peak_times(cnum) = time_closest_pk;
            
        end % end loop through contacts on one probe
        raw_signal = raw_signal(subwin_st:subwin_end,:);
        zscore_signal = zscore_signal(subwin_st:subwin_end,:);
        filt_signal = filt_signal(subwin_st:subwin_end,:);

        % % plot filt_signal and display it
        colors = [
            0.76 0.86 0.98;
            %0.72 0.83 0.96;
            %0.68 0.79 0.94;
            %0.64 0.76 0.92;
            0.60 0.72 0.90;
            %0.56 0.69 0.88;
            %0.52 0.65 0.86;
            %0.48 0.62 0.84;
            0.44 0.58 0.82;
            %0.40 0.55 0.80;
            %0.36 0.51 0.78;
            %0.32 0.48 0.76;
            0.28 0.44 0.74;
            %0.24 0.41 0.72;
            %0.20 0.37 0.70;
            %0.16 0.34 0.68;
            0.12 0.30 0.66;
            0.06 0.26 0.63];
        % colors = [
        %     0 0 0; % TEMP - delete this row after 1.25.26
        %     0.76 0.86 0.98;
        %     0.72 0.83 0.96;
        %     0.68 0.79 0.94;
        %     0.64 0.76 0.92;
        %     0.60 0.72 0.90;
        %     0.56 0.69 0.88;
        %     % 0.52 0.65 0.86;
        %     % 0.48 0.62 0.84;
        %     % 0.44 0.58 0.82;
        %     0.40 0.55 0.80;
        %     % 0.36 0.51 0.78;
        %     % 0.32 0.48 0.76;
        %     % 0.28 0.44 0.74;
        %     0.24 0.41 0.72;
        %     % 0.20 0.37 0.70;
        %     % 0.16 0.34 0.68;
        %     % 0.12 0.30 0.66;
        %     0.06 0.26 0.63];
        % fg1 = figure; 
        % set( fg1, 'position', [ 88  100  1400  600 ] ) %[ 88  1593  1250  420 ]
        % for cnum = 1:numel(cur_elec_contact_ind)
        %     plot(filt_signal(:,cnum),'linewidth',2,'color',colors(cnum,:)) % -1024:-1
        %     hold on
        %     xlim([subwin_st,subwin_end]); 
        %     %ylim([-0.5, 0.5]);
        %     xlabel('samples from alignment'); 
        %     ylabel('zscored amplitude');
        %     yl1 = ylim;
        %     title(sprintf("%s Elec %s (%d-%d)",alignment,cur_letter, subwin_st, subwin_end))
        % end
        % 
        % close;

        % check for traveling waves in peak_times
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
                % check whether traveling wave is in same direction,
                % involving same contacts (before 1.24.26: ANY overlap in contacts is acceptable)

                ct_incr = [idx-increase_count idx-1];
                % trial_avg_direction = traveling_ERP_table{row,"Direction"};
                % trial_avg_cts = traveling_ERP_table{row,"Contacts"};
                % overlap = (ct_incr(1) >= trial_avg_cts(1) && ct_incr(1) <= trial_avg_cts(2)) || (ct_incr(2) >= trial_avg_cts(1) && ct_incr(2) <= trial_avg_cts(2));
                

                % 1.24.26 -- require that all (4+) contacts
                % participating in traveling wave in trial-avg
                % are invovled in the single-trial wave
                    % HARD-CODED at top of section
                overlap = (ct_incr(1) <= trial_avg_cts(1)) && (ct_incr(2) >= trial_avg_cts(2));
                if (trial_avg_direction == 1) & (overlap)
                    single_trial_count = single_trial_count + 1;
                end
            end
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

                ct_decr = [idx-decrease_count idx-1];
                % trial_avg_direction = traveling_ERP_table{row,"Direction"};
                % trial_avg_cts = traveling_ERP_table{row,"Contacts"};
                % overlap = (ct_decr(1) >= trial_avg_cts(1) && ct_decr(1) <= trial_avg_cts(2)) || (ct_decr(2) >= trial_avg_cts(1) && ct_decr(2) <= trial_avg_cts(2));
                overlap = (trial_avg_cts(1) >= ct_decr(1)) && (trial_avg_cts(2) <= ct_decr(2));
                if (trial_avg_direction == -1) & (overlap)
                    single_trial_count = single_trial_count + 1;
                    single_trials = [single_trials; iter];


                    % PLOT
                    fg1 = figure; 
                    set( fg1, 'position', [ 88  100  1100  700 ] ) %[ 88  1593  1250  420 ]
                    tiledlayout(2,1,'TileSpacing','compact');

                    nexttile;
                    for cnum = trial_avg_cts(1):trial_avg_cts(2)
                        plot(filt_signal(:,cnum),'linewidth',2,'color',colors(cnum-trial_avg_cts(1)+1,:)) % -1024:-1
                        hold on
                        xlim([0,200]);
                        %xlim([subwin_st,subwin_end]); 
                        %ylim([-0.5, 0.5]);
                        xlabel('samples from alignment'); 
                        ylabel('zscored amplitude');
                        yl1 = ylim;
                        title(sprintf("%s Elec %s (%d-%d)",alignment,cur_letter, subwin_st, subwin_end))
                    end

                    nexttile;
                    for cnum = trial_avg_cts(1):trial_avg_cts(2)
                        plot(zscore_signal(:,cnum),'linewidth',2,'color',colors(cnum-trial_avg_cts(1)+1,:)) % -1024:-1
                        hold on
                        xlim([0,200]);
                        %xlim([subwin_st,subwin_end]); 
                        %ylim([-0.5, 0.5]);
                        xlabel('samples from alignment'); 
                        ylabel('zscored amplitude');
                        yl1 = ylim;
                        title("Z-SCORED (NO FILT)")
                    end

                    % nexttile;
                    % for cnum = 1:numel(cur_elec_contact_ind)
                    %     plot(raw_signal(:,cnum),'linewidth',2,'color',colors(cnum,:)) % -1024:-1
                    %     hold on
                    %     xlim([0,200]);
                    %     %xlim([subwin_st,subwin_end]); 
                    %     %ylim([-0.5, 0.5]);
                    %     xlabel('samples from alignment'); 
                    %     ylabel('amplitude');
                    %     yl1 = ylim;
                    %     title("RAW")
                    % end
                    % 
                    % %1.25.26
                    fname = sprintf("%s/waveform_plots/3.3.26/trial%d_match_trialavg.jpg",out_dir,iter);
                    print(fg1,'-djpeg',fname)
                    % 
                    close;

                end
            end
            decrease_count = 1; % reset
        
            while idx <= numel(peak_times) && peak_times(idx) == prev_latency
                prev_latency = peak_times(idx);
                idx = idx + 1;
            end
        
        end

    end % end loop through events

end % end loop through sessions

fprintf("Elec N (fourth inspection): %d single trials with waves that match main one in trial-avg \n", single_trial_count)
% END RUN HERE 1/24/26

%% 2.22.26 - Visualize a few single trials that DON'T have matching traveling waves
% to compare ground vs neighbor-average -- how the latter affects phase organization
% Elec N (trials/events NOT IN 'single_trials' array)
% Run this section right after above section

for iter = 1:numel(align_times)  
    rand_event = iter;

    fg1 = figure; 
    set( fg1, 'position', [ 88  100  1100  700 ] ) %[ 88  1593  1250  420 ]
    tiledlayout(2,1,'TileSpacing','compact');

    for refnum = 1:numel(references)
        reference = references{refnum};
        ref_name = ref_names{refnum};

        raw_signal = []; zscore_signal = [];
        filt_signal = []; % each column = signal of one contact
        peak_times = zeros(numel(cur_elec_contact_ind),1);
        for cnum = 1:numel(cur_elec_contact_ind)
            cur_name = cur_elec_contact_names(cnum);
            contact = cur_elec_contact_ind(cnum);
            
            % load real-valued signal (separate_channel_files -> filter)
            if strcmp(reference,'Ground')
                load(sprintf('%s/%s/separate_channel_files/%s_MAD_SES%d_ch%03d.mat',data_base_dir,subject_ID,subject_ID,sesnum,contact)); 
            elseif strcmp(reference,'neighbor_average')
                load(sprintf('%s/%s/separate_channel_files/%s/%s_MAD_SES%d_ch%03d.mat',data_base_dir,subject_ID,reference,subject_ID,sesnum,contact));
            end
            
            raw = data; raw_event = raw(align_times(rand_event)-win(1):align_times(rand_event)+win(2)-1);
            zscored_event = (raw_event - mean(raw_event)) ./ std(raw_event);
            [b,a] = butter( 4, [5 40] ./ (Fs/2) ); data = filtfilt( b, a, data );
            data_event = data(align_times(rand_event)-win(1):align_times(rand_event)+win(2)-1);
            data_event = ( data_event - mean(data_event) ) ./ std(data_event);
            
            raw_signal = [raw_signal raw_event];
            zscore_signal = [zscore_signal zscored_event];
            filt_signal = [filt_signal data_event];
    
            cur_ct_pk_times = [];
            if peak_wave
                cur_ct_pk_times = find(sign(diff([data_event; data_event(end)]))<0&sign(diff([data_event(1); data_event]))>0);
            else
                cur_ct_pk_times = find(sign(diff([data_event; data_event(end)]))>0&sign(diff([data_event(1); data_event]))<0);
            end
    
            [~,idx_mindiff] = min(abs(cur_ct_pk_times - time_max_pk));
            time_closest_pk = cur_ct_pk_times(idx_mindiff);
            
            peak_times(cnum) = time_closest_pk;
            
        end % end loop through contacts on one probe
        raw_signal = raw_signal(subwin_st:subwin_end,:);
        zscore_signal = zscore_signal(subwin_st:subwin_end,:);
        filt_signal = filt_signal(subwin_st:subwin_end,:);
    
        nexttile; 
        for cnum = trial_avg_cts(1):trial_avg_cts(2)
            plot(filt_signal(:,cnum),'linewidth',2,'color',colors(cnum-trial_avg_cts(1)+1,:)) % -1024:-1
            hold on
            xlim([0,200]);
            %xlim([subwin_st,subwin_end]); 
            %ylim([-0.5, 0.5]);
            xlabel('samples from alignment'); 
            ylabel('zscored amplitude');
            yl1 = ylim;
            title(sprintf("%s %s Elec %s (%d-%d) Trial %d",ref_name, align_name,cur_letter, subwin_st, subwin_end, iter))
        end
    
    end % end loop through references

    fname = sprintf("%s/waveform_plots/2.21.26/trial%d_match_trialavg.jpg",out_dir,iter);
    %print(fg1,'-djpeg',fname)
    close;

end % end loop through events

%% 1.29.26 -- Detect and count traveling waves in aligned trial-avg and jittered trial-avg


%% What percent of lists of 18 integers 1-100 contain contiguous runs that are increasing/decreasing?
