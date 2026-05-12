% ERP -> traveling wave detection

addpath('/media/Data/Human_Intracranial_MAD/_toolbox/spectral-analysis');
subject_IDs = {'EMU001','EMU024', 'EMU025','EMU030','EMU037','EMU038','EMU039','EMU041','EMU047','EMU051'};
reference = 'Ground'; 
datapath = '/media/Data/Human_Intracranial_MAD/';
data_base_dir = sprintf('%s1_formatted',datapath);
addpath('/media/Data/Human_Intracranial_MAD/analysis/PAC_code/CF_Coupling/generalized phase')

%% 1. Save list of p-values for each peak and trough
% Compare aligned trial-avg signal to time-jittered
alignment = 'inspection';

for subject_num = 1 %1:numel(subject_IDs)
    
    subject_ID = subject_IDs{subject_num};
    
    if strcmp(subject_ID,'EMU001')
        Fs = 1000;
        win = [0 500];
    else
        Fs = 2048;
        win = [0 1024];
    end

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

        j = 1;
        if ~contains(cur_letter, "'")
            while ~contains(elec_name{j},cur_letter) || contains(elec_name{j},"'"), j = j+1; end
        else
            while ~contains(elec_name{j},cur_letter), j = j+1; end
        end
        cur_name = elec_name{j};
        cur_elec_contact_names = []; % stores name of contacts on current electrode
        cur_elec_contact_ind = []; % stores indices of contacts on curr elec
        while strcmp(cur_name(1:2),cur_letter) || (strcmp(cur_name(1),cur_letter) && ~contains(cur_name, "'"))
            cur_elec_contact_names = [cur_elec_contact_names; convertCharsToStrings(cur_name)];
            cur_elec_contact_ind = [cur_elec_contact_ind; j];
            
            j = j+1; 
            if j > numel(elec_name), break; end
            cur_name = elec_name{j};
        end

        contacts = cur_elec_contact_ind;

        %alignments = ["inspection"];
        %alignments = ["outcome"];
        %alignments = ["probability","amount"];
        %alignments = ["first_unique_attribute", "fourth_unique_attribute"];
        %alignments = ["second_unique_attribute", "third_unique_attribute"];
        %alignments = ["likely_outcome","unlikely_outcome"];
        %alignments = ["positive", "negative"];
        alignments = ["inspection","single_opt_first_inspection","full_single_opt_info"];

        for alignment = alignments

            out_dir = sprintf('/media/Data/Human_Intracranial_MAD/analysis/Cross_Freq_Coupling/%s/%s/synchrony_data/generalized_phase',subject_ID, reference);  
            % final_out_dir = sprintf('%s/waveform_plots/trial_averaged/%s',out_dir, alignment);
            % if ~exist(final_out_dir, 'dir'), mkdir(final_out_dir); end
        
            % data_out_dir = sprintf('%s/generalized_phase/data_filtered_trial_avg/%s',out_dir, alignment);
            % if ~exist(data_out_dir, 'dir'), mkdir(data_out_dir); end

            for cnum = 1:numel(contacts) 
                contact = contacts(cnum);
                cur_name = cur_elec_contact_names(cnum);
                %channel = cur_elec_contact_ind(cnum);
                
                num_events = 0;
                sum_events = zeros(max(win),1);
                %sum_squared_events = zeros(max(win),1);
                for sesnum = 1:num_sessions
                    load(sprintf('%s/%s/%s_MAD_SES%d_Setup.mat',data_base_dir,subject_ID,subject_ID,sesnum));
                    
                    [align_times,trial_numbers] = get_align_times(filters, trial_times, trial_words, alignment);
    
                    remove_ind = isnan(align_times);
                    align_times(isnan(align_times)) = [];
                    trial_numbers(remove_ind) = [];
                    
                    align_times = round(align_times*Fs);  

                    if strcmp(subject_ID,'EMU001')
                        [align_times_inspect,trial_numbers_inspect] = get_align_times(filters, trial_times, trial_words, '2opt_2att_inspection');

                        remove_ind = isnan(align_times_inspect);
                        align_times_inspect(isnan(align_times_inspect)) = [];
                        trial_numbers_inspect(remove_ind) = [];

                        align_times_inspect = round(align_times_inspect*Fs); 

                        align_times = intersect(align_times, align_times_inspect);
                    end
                    
                    load(sprintf('%s/%s/separate_channel_files/%s_MAD_SES%d_ch%03d.mat',data_base_dir,subject_ID,subject_ID,sesnum,contact)); % for saving phase for ALL contacts
                    
                    % broadband filter
                    [b,a] = butter( 4, [5 40] ./ (Fs/2) ); filtered = filtfilt( b, a, data ); 
                    data = filtered; 

                    for ii = 1:numel(align_times)
        
                        data_event = data((round(align_times(ii)-win(1)):round(align_times(ii)+win(2)-1)));
                        data_event = ( data_event - mean(data_event) ) ./ std(data_event);
                        sum_events = sum_events + data_event;      

                    end % end loop through events
                    num_events = num_events + numel(align_times);
                end % end loop through sessions

                avg_events_aligned = sum_events./ num_events;

                % FIND ALL PEAKS AND TROUGHS
                % Construct complex signal
                lp = 5; xgp_avg = generalized_phase_vector( avg_events_aligned, Fs, lp );  
                % Get angle
                angle_avg = angle(xgp_avg);
                angle_posneg = angle_avg;
                angle_posneg(angle_posneg >= 0) = 1;
                angle_posneg(angle_posneg < 0) = -1;

                peaks_idx_align = find(diff(angle_posneg) > 0);
                troughs_idx_align = find(diff(angle_posneg) < 0);
                peaks_align = avg_events_aligned(peaks_idx_align);
                troughs_align = avg_events_aligned(troughs_idx_align);

                iterations = 1000;
                peak_trough_distr = [];

                for iter = 1:iterations

                    sum_events = zeros(max(win),1);
                    num_events = 0;
                    for sesnum = 1:num_sessions
                        load(sprintf('%s/%s/%s_MAD_SES%d_Setup.mat',data_base_dir,subject_ID,subject_ID,sesnum));
    
                        [align_times,trial_numbers] = get_align_times(filters, trial_times, trial_words, alignment);
    
                        remove_ind = isnan(align_times);
                        align_times(isnan(align_times)) = [];
                        trial_numbers(remove_ind) = [];
    
                        align_times = round(align_times*Fs);

                        if strcmp(subject_ID,'EMU001')
                            [align_times_inspect,trial_numbers_inspect] = get_align_times(filters, trial_times, trial_words, '2opt_2att_inspection');
    
                            remove_ind = isnan(align_times_inspect);
                            align_times_inspect(isnan(align_times_inspect)) = [];
                            trial_numbers_inspect(remove_ind) = [];
    
                            align_times_inspect = round(align_times_inspect*Fs); 
    
                            align_times = intersect(align_times, align_times_inspect);
                        end
    
                        load(sprintf('%s/%s/separate_channel_files/%s_MAD_SES%d_ch%03d.mat',data_base_dir,subject_ID,subject_ID,sesnum,contact)); % for saving phase for ALL contacts
    
                        % RNG for random time offsets
                        rng(iter+sesnum, 'twister') % seed = iter+sesnum here
                        rand_nums = rand(numel(align_times),1);
    
                        % broadband filter
                        [b,a] = butter( 4, [5 40] ./ (Fs/2) ); filtered = filtfilt( b, a, data ); 
                        data = filtered;
    
                        %sum_events = zeros(max(win),1);
                        for ii = 1:numel(align_times)
    
                            offset_sec = -0.25 + 0.5 * rand_nums(ii);
                            offset_samples = offset_sec * Fs;
    
                            data_event = data((round(align_times(ii)-win(1)+offset_samples):round(align_times(ii)+win(2)-1+offset_samples)));
                            data_event = (data_event - mean(data_event))./ std(data_event); % zscore - TEMP
                            sum_events = sum_events + data_event;    
    
                        end % end loop through events
                        num_events = num_events + numel(align_times);
                    end % end loop through sesnums
                    avg_events_jittered = sum_events ./ num_events;

                    % % FIND PEAKS AND TROUGHS
                    lp = 5; xgp_avg_jit = generalized_phase_vector( avg_events_jittered, Fs, lp );  
                    angle_avg_jit = angle(xgp_avg_jit);
                    angle_posneg = angle_avg_jit;
                    angle_posneg(angle_posneg >= 0) = 1;
                    angle_posneg(angle_posneg < 0) = -1;

                    peaks_idx_jit = find(diff(angle_posneg) > 0);
                    troughs_idx_jit = find(diff(angle_posneg) < 0);
                    peaks_jit = avg_events_jittered(peaks_idx_jit);
                    troughs_jit = avg_events_jittered(troughs_idx_jit);
    
                    % add peak and tough amplitudes to respective distributions
                    peak_trough_distr = [peak_trough_distr; peaks_jit; troughs_jit];

                end % end loop through iterations

                peak_trough_distr = sort(peak_trough_distr,'ascend');

                pv_table = cell(0,4); % one table for each contact
                for peak_num = 1:numel(peaks_align)
                    peak = peaks_align(peak_num); % zscored peak
                    
                    peak_idx = peaks_idx_align(peak_num);
                    % find the trough right before and trough right after
                    idx_trough_before = troughs_idx_align(troughs_idx_align < peak_idx);
                    if ~isempty(idx_trough_before)
                        idx_trough_before = idx_trough_before(end);
                    else
                        idx_trough_before = -1;
                    end
                    idx_trough_after = troughs_idx_align(troughs_idx_align > peak_idx);
                    if ~isempty(idx_trough_after)
                        idx_trough_after = idx_trough_after(1);
                    else
                        idx_trough_after = -1;
                    end
                    % FIND P-VALUE
                    [~,pv_idx] = min(abs(peak_trough_distr-peak));
                    pv = 1 - (pv_idx/numel(peak_trough_distr));
                    pv_table(end+1,:) = {"peak", peak_idx, peak, pv};
                end
                    
                for trough_num = 1:numel(troughs_align)
                    trough = troughs_align(trough_num);
                    
                    trough_idx = troughs_idx_align(trough_num);
                    idx_peak_before = peaks_idx_align(peaks_idx_align < trough_idx);
                    if ~isempty(idx_peak_before)
                        idx_peak_before = idx_peak_before(end);
                    else
                        idx_peak_before = -1;
                    end
                    
                    idx_peak_after = peaks_idx_align(peaks_idx_align > trough_idx);
                    if ~isempty(idx_peak_after)
                        idx_peak_after = idx_peak_after(1);
                    else
                        idx_peak_after = -1;
                    end
                    % FIND P-VALUE
                    [~,pv_idx] = min(abs(peak_trough_distr-trough));
                    pv = pv_idx/numel(peak_trough_distr);
                    pv_table(end+1,:) = {"trough", trough_idx, trough, pv};
                end
                pv_table = cell2table(pv_table, "VariableNames",["Peak/Trough","Idx","Amplitude","p-value"]);
                if ~exist(sprintf('%s/ERP_detection/%s',out_dir,alignment), 'dir'), mkdir(sprintf('%s/ERP_detection/%s',out_dir,alignment)); end
                fname_table = sprintf('%s/ERP_detection/%s/pvals_%s_%03d.mat', out_dir,alignment,cur_name, contact); 
                save(fname_table, "pv_table")

            end % end loop through contacts
        
        end % end loop through alignments
        disp("elec %s tables saved", cur_letter)

    end % end loop through electrodes

end %end loops through subj

%% 2. Create ERP_sig_table
% Load p-values of all peaks and troughs to test for significance and apply fdr correction

alignment = 'inspection';

for subject_num = 1 %[1 5 7 8] 
    
    subject_ID = subject_IDs{subject_num};
    
    if strcmp(subject_ID,'EMU001')
        Fs = 1000;
        win = [0 500];
    else
        Fs = 2048;
        win = [0 1024];
    end

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

    % get list of elec/probe letters
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

        j = 1;
        if ~contains(cur_letter, "'")
            while ~contains(elec_name{j},cur_letter) || contains(elec_name{j},"'"), j = j+1; end
        else
            while ~contains(elec_name{j},cur_letter), j = j+1; end
        end
        cur_name = elec_name{j};
        cur_elec_contact_names = []; % stores name of contacts on current electrode
        cur_elec_contact_ind = []; % stores indices of contacts on curr elec
        while strcmp(cur_name(1:2),cur_letter) || (strcmp(cur_name(1),cur_letter) && ~contains(cur_name, "'"))
            cur_elec_contact_names = [cur_elec_contact_names; convertCharsToStrings(cur_name)];
            cur_elec_contact_ind = [cur_elec_contact_ind; j];
            
            j = j+1; 
            if j > numel(elec_name), break; end
            cur_name = elec_name{j};
        end

        contacts = cur_elec_contact_ind;

        ERP_table = table(); 

        % create a table of significant ERPs on each contact
        alignments = ["inspection","single_opt_first_inspection","full_single_opt_info"];

        for alignment = alignments
            out_dir = sprintf('/media/Data/Human_Intracranial_MAD/analysis/Cross_Freq_Coupling/%s/%s/synchrony_data/generalized_phase',subject_ID, reference);  
            
            fprintf("%s #sig ERPs Elec %s\n", alignment, cur_letter)
            num_sig = 0;
            for cnum = 1:numel(contacts) 
                contact = contacts(cnum);
                cur_name = cur_elec_contact_names(cnum);

                % load saved p-values
                fname_table = sprintf('%s/ERP_detection/%s/pvals_%s_%03d.mat', out_dir,alignment,cur_name, contact); 
                load(fname_table) % pv_table is loaded
                pvals = pv_table(:,"p-value");
                pvals = table2array(pvals);

                % run fdr_hb
                q = 0.05; method = 'dep'; report = 'no';
                [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(pvals,q,method,report);

                % we want to save h and adj_p
                fname_pvs = sprintf('%s/ERP_detection/%s/sig_ERPs_%s_%03d.mat', out_dir,alignment,cur_name, contact); 
                sig_pv = addvars(pv_table, h, 'NewVariableNames', 'Sig');
                sig_pv = addvars(sig_pv, adj_p, 'NewVariableNames', 'Adj p-val');
                % save(fname_pvs, "sig_pv")

                fprintf("%s: %d\n", cur_name, nnz(h))
                num_sig = num_sig + nnz(h);

                cur_elec_col = repmat(sprintf("%s",cur_name), height(sig_pv), 1);
                cur_ct_sig_pv = addvars(sig_pv, cur_elec_col, 'NewVariableNames','Contact','Before', 'Peak/Trough');
                cur_align_col = repmat(sprintf("%s",alignment),height(sig_pv),1);
                cur_ct_sig_pv = addvars(cur_ct_sig_pv, cur_align_col,'NewVariableNames','Alignment','Before', 'Peak/Trough');
                ERP_table = [ERP_table; cur_ct_sig_pv];

            end % end loop through contacts
            fprintf("Total sig ERPs %s: %d\n", cur_letter, num_sig)
        end % end loop through alignments

        % save ERP_table
        fname_table = sprintf('%s/fdr_ERP_sig_table_%s.mat', out_dir,cur_letter); 
        save(fname_table, "ERP_table")

    end % end loop through electrodes
end % end loop through subj

%% 3. Iterate through ERP_table to find traveling waves surrounding each ERP

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

            % Create new table for saving traveling waves detected
            traveling_ERP_table = cell(0,10);

            for alignment = alignments

                per_contact_wave_counts = zeros(numel(cur_elec_contact_ind),1);
                per_contact_ERP_counts = zeros(numel(cur_elec_contact_ind),1);

                travel_ERP_count = 0;
                
                out_dir = sprintf('/media/Data/Human_Intracranial_MAD/analysis/Cross_Freq_Coupling/%s/%s/synchrony_data/generalized_phase',subject_ID, reference);  
                final_out_dir = sprintf('%s/waveform_plots/trial_averaged/%s',out_dir, alignment);
                data_out_dir = sprintf('%s/data_filtered_trial_avg/%s',out_dir, alignment);
            
                ERP_peaks_align = ERP_peaks(strcmp(ERP_peaks{:,"Alignment"},alignment),:);
                ERP_troughs_align = ERP_troughs(strcmp(ERP_troughs{:,"Alignment"},alignment),:);

                % Delete all rows in ERP_peaks_align and ERP_troughs_align that are not significant
                ERP_peaks_align(ERP_peaks_align{:,"Sig"}==0,:) = [];
                ERP_troughs_align(ERP_troughs_align{:,"Sig"}==0,:) = [];

                % 10/30/25 
                % (boolean) true = this row's ERP participates in traveling wave; 
                %           false = does not participate in traveling wave
                trav_pks_bool = zeros(height(ERP_peaks_align),1);
                trav_tghs_bool = zeros(height(ERP_troughs_align),1);

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
                                p = polyfit(1:increase_count,peak_times(idx-increase_count:idx-1),1);
                                incr_speed = 3.5 / (p(1) * (1000/2048));
                                circlin_corr = zeros(height(angles_subwin),2);
                                for t = 1:height(angles_subwin)
                                    [cc, pv] = circ_corrcl(1:increase_count,angles_subwin(t,idx-increase_count:idx-1));
                                    circlin_corr(t,1) = cc; circlin_corr(t,2) = pv;
                                end
                                traveling_ERP_table(end + 1, :) = {alignment, "peak", subwin_start, subwin_end, [idx-increase_count idx-1], 1, peak_times, incr_speed, circlin_corr, ERP_table{row,1}};
                                
                                % Check whether detected wave encompasses the contact 
                                % where this row's ERP was detected
                                if (cnum >= idx-increase_count) & (cnum <= idx-1) % added 10/10/25
                                    travel_ERP_count = travel_ERP_count + 1;
                                    per_contact_wave_counts(cnum) = per_contact_wave_counts(cnum) + 1;
                                    trav_pks_bool(row) = 1;
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

                                traveling_ERP_table(end + 1, :) = {alignment, "peak", subwin_start, subwin_end, [idx-decrease_count idx-1], -1, peak_times, decr_speed, circlin_corr, ERP_table{row,1}};
                                
                                % Check whether detected wave encompasses the contact 
                                % where this row's ERP was detected
                                if (cnum >= idx-decrease_count) & (cnum <= idx-1) % added 10/10/25
                                    travel_ERP_count = travel_ERP_count + 1;
                                    per_contact_wave_counts(cnum) = per_contact_wave_counts(cnum) + 1;
                                    trav_pks_bool(row) = 1;
                                end
                            end

                            decrease_count = 1; % reset
                        
                            while idx <= numel(peak_times) && peak_times(idx) == prev_latency
                                prev_latency = peak_times(idx);
                                idx = idx + 1;
                            end
                        
                        end

                    end % end loop through rows in ERP_peaks_align
                
                end % end IF peaks non-empty

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

                            decrease_count = 1; % reset
                        
                            while idx <= numel(peak_times) && peak_times(idx) == prev_latency
                                prev_latency = peak_times(idx);
                                idx = idx + 1;
                            end
                        
                        end
    
                    end % end loop through rows of trough ERP table
                    
                end % end IF trough table non-empty

            end % end loop through alignments

            out_dir = sprintf('/media/Data/Human_Intracranial_MAD/analysis/Cross_Freq_Coupling/%s/%s/synchrony_data/generalized_phase',subject_ID, reference);
            table_out_dir = sprintf('%s/phase_progression_detection/trialavg',out_dir);
            if ~exist(table_out_dir, "dir"), mkdir(table_out_dir); end
            traveling_ERP_table = cell2table(traveling_ERP_table,"VariableNames",["Alignment","Peak/Trough","Subwin St","Subwin End","Contacts","Direction","Latencies","Speed (m/s)","Corrs", "Ref Ct"]);
            fname_table = sprintf('%s/%s_trialavg_travERP_table.mat', table_out_dir, cur_letter); 
            save(fname_table, "traveling_ERP_table")

        end % end loop through elecs

    %end % end loop through alignments
    
end % end loop through subjects