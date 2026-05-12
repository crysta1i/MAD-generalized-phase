%% 0. Split white and gray matter

% for subj

    % for elec

        wm_gm_labels = strings(numel(cur_elec_contact_ind),1);
        for cnum = 1:numel(cur_elec_contact_ind)
            if strcmpi(elec_area(cur_elec_contact_ind(cnum)),"WM")
                wm_gm_labels(cnum) = "WM";
            elseif strcmpi(elec_area(cur_elec_contact_ind(cnum)),"undefined") || strcmpi(elec_area(cur_elec_contact_ind(cnum)),"OUT")
                wm_gm_labels(cnum) = "other";
            elseif strcmpi(elec_area(cur_elec_contact_ind(cnum)),"bolt")
                wm_gm_labels(cnum) = "other";
            else
                wm_gm_labels(cnum) = "GM";
            end
        end

        cur_wm_gm_chunks = [];
        cur_wm_gm_chunk_names = [];

        cnum = 2;
        while cnum <= numel(cur_elec_contact_names)
            contact = cur_elec_contact_ind(cnum);
            prev_contact = cur_elec_contact_ind(cnum-1);
            cur_name = cur_elec_contact_names(cnum);
            cur_area = wm_gm_labels(cnum);
            prev_area = wm_gm_labels(cnum-1);

            if ~strcmpi(cur_area,prev_area)
                cur_wm_gm_chunks = [cur_wm_gm_chunks; contact];
                cur_wm_gm_chunk_names = [cur_wm_gm_chunk_names; cur_name];
            end

            cnum = cnum + 1;
        end

%% 1. Distribution of inter-peak times

% For contacts on one electrode (1d waves)


% TODO: can build this code on top of detect_peak_timing code...
    % Section 3 of detect_peak_timing: not constrained to detected ERPs
    % Section 4 of detect_peak_timing: constrained to detected ERPs

% Separate wm from gm
% -- add speed (or inter-peak/inter-trough time) to distribution



% Need own script for this part (here): check for trav waves within chunks
all_probe_speed_distr = [];
all_probe_gm_speed_distr = [];
all_probe_wm_speed_distr = [];

% for each electrode

    % load trial-avg data (real and xgp) for each contact, concatenate
    % create peak_times / trough_times array
    
    gap_distr_wm = [];
    gap_distr_gm = [];
    gap_distr_all = []; % includes gaps btwn peaks of a contact in gm and a contact in wm

    travel_ERP_count = 0;
    trav_wave_idx = []; 
    % a N x 2 matrix where N = number of detected traveling waves
    % 1st column: index (in cur_elec_contact_names) where traveling wave STARTS
    % 2nd column: index where traveling wave ENDS (inclusive)
    % Does not distinguish direction of wave -- all are added to this array

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
            trav_wave_idx = [trav_wave_idx; [idx-increase_count idx-1]];
            travel_ERP_count = travel_ERP_count + 1;
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
            trav_wave_idx = [trav_wave_idx; [idx-decrease_count idx-1]];
            travel_ERP_count = travel_ERP_count + 1;
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
    
    

% end % end loop through electrodes

gm_speed_distr = [];
wm_speed_distr = [];

%% 2. Distribution of best-fit line slopes
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
    wm_gm_elec_area = strings(numel(elec_area),1);
    for c = 1:numel(elec_area)
        if strcmpi(elec_area(c),"WM")
            wm_gm_elec_area(c) = "WM";
        elseif strcmpi(elec_area(c),"undefined") || strcmpi(elec_area(c),"OUT")
            wm_gm_elec_area(c) = "other";
        elseif strcmpi(elec_area(c),"bolt")
            wm_gm_elec_area(c) = "other";
        else
            wm_gm_elec_area(c) = "GM";
        end
    end


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

    distr_speeds = []; % all electrodes, all three alignments
    distr_speeds_wm = []; distr_speeds_gm = [];
    distr_speeds_inspect1 = []; distr_speeds_inspect2 = []; distr_speeds_inspect = [];
    amplitudes = []; amplitudes_gm = [];
    % scatterplot (1) distr_speeds vs amplitudes, (2) distr_speeds_gm vs amplitudes_gm
    for e = 1:numel(elec_letters)
        cur_letter = elec_letters(e);

        % get cur_elec_contact_ind/names
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

        % load traveling_ERP_table
        fname_table = sprintf('%s/%s_trialavg_travERP_table.mat', table_out_dir, cur_letter); 
        load(fname_table, "traveling_ERP_table")

        for row = 1:height(traveling_ERP_table)
            direction = traveling_ERP_table{row, "Direction"};
            alignment = traveling_ERP_table{row, "Alignment"};
            contacts = traveling_ERP_table{row, "Contacts"};
            ct_st = cur_elec_contact_ind(contacts(1)); 
            ct_end = cur_elec_contact_ind(contacts(2));
            
            amplitudes = [amplitudes; traveling_ERP_table{row, 10}];
            distr_speeds = [distr_speeds; traveling_ERP_table{row,8}];
            if all(wm_gm_elec_area(ct_st:ct_end) == "WM")
                distr_speeds_wm = [distr_speeds_wm; traveling_ERP_table{row, 8}];
            elseif all(wm_gm_elec_area(ct_st:ct_end) == "GM")
                distr_speeds_gm = [distr_speeds_gm; traveling_ERP_table{row, 8}];
                amplitudes_gm = [amplitudes_gm; traveling_ERP_table{row, 10}];
            end

            if all(wm_gm_elec_area(ct_st:ct_end) == "GM") && strcmp(alignment,"single_opt_first_inspection")
                distr_speeds_inspect1 = [distr_speeds_inspect1; traveling_ERP_table{row, 8}];
            elseif all(wm_gm_elec_area(ct_st:ct_end) == "GM") && strcmp(alignment,"full_single_opt_info")
                distr_speeds_inspect2 = [distr_speeds_inspect2; traveling_ERP_table{row, 8}];
            elseif all(wm_gm_elec_area(ct_st:ct_end) == "GM") && strcmp(alignment,"inspection")
                distr_speeds_inspect = [distr_speeds_inspect; traveling_ERP_table{row, 8}];
            end
            
            % if strcmp(alignment,"single_opt_first_inspection")
            %     distr_speeds_inspect1 = [distr_speeds_inspect1; traveling_ERP_table{row, 8}];
            % elseif strcmp(alignment,"full_single_opt_info")
            %     distr_speeds_inspect2 = [distr_speeds_inspect2; traveling_ERP_table{row, 8}];
            % else
            %     distr_speeds_inspect = [distr_speeds_inspect; traveling_ERP_table{row, 8}];
            % end
        end % end loop through rows of trav ERP table

    end % end loop through electrodes

    % figure;
    % tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    % nexttile;
    % plot(amplitudes_gm, abs(distr_speeds_gm),'.')
    % title("Velocity vs Max Abs Amplitude - GM only")
    % nexttile;
    % plot(amplitudes, abs(distr_speeds),'.')
    % title("Velocity vs Max Abs Amplitude")

    % figure;
    % tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
    % 
    % nexttile;
    % histogram(distr_speeds_inspect,'BinWidth', 0.1)
    % title("Inspection")
    % xlim([-2,2])
    % 
    % nexttile;
    % histogram(distr_speeds_inspect1,'BinWidth', 0.1)
    % title("1st inspection of opt")
    % xlim([-2,2])
    % 
    % nexttile;
    % histogram(distr_speeds_inspect2,'BinWidth', 0.1)
    % title("2nd inspection of opt")
    % xlim([-2,2])
    
    figure;
    histogram(distr_speeds(abs(distr_speeds)<=3),'BinWidth', 0.1)
    title(sprintf("Distribution of Traveling ERP Velocities (%s)",subject_ID))
    xlim([-3,3])
    xlabel("velocity (m/s)")
    ylabel("count")

    fprintf("Mean velocity: %f \n", mean(distr_speeds(abs(distr_speeds)<=3)))
    fprintf("Median velocity: %f \n", median(distr_speeds(abs(distr_speeds)<=3)))
    fprintf("Mean speed: %f \n", mean(abs(distr_speeds(abs(distr_speeds)<=3))))
    fprintf("Median speed: %f \n", median(abs(distr_speeds(abs(distr_speeds)<=3))))

    % nexttile;
    % histogram(distr_speeds,'BinWidth', 0.1)
    % title("Distr of traveling ERP speeds (m/s)")
    % xlim([-2,2])
    % 
    % nexttile;
    % histogram(distr_speeds_gm, 'BinWidth',0.1)
    % title("Gray Matter")
    % xlim([-2,2])

    % nexttile;
    % histogram(distr_speeds_wm, 'BinWidth',0.1)
    % title("White Matter")
    % xlim([-2,2])

end % end loop through subj

%% TEMP: RUN THIS ONCE FOR EACH SUBJECT
% Add column of max_amp to the traveling_ERP_table
% DONE: EMU038

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

    distr_speeds = []; % all electrodes, all three alignments
    distr_speeds_wm = []; distr_speeds_gm = [];
    distr_speeds_inspect1 = []; distr_speeds_inspect2 = []; distr_speeds_inspect = [];
    for e = 1:numel(elec_letters)
        cur_letter = elec_letters(e);

        % get cur_elec_contact_ind/names
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

        % load traveling_ERP_table
        fname_table = sprintf('%s/%s_trialavg_travERP_table.mat', table_out_dir, cur_letter); 
        load(fname_table, "traveling_ERP_table")

        max_amps = []; % column
        for row = 1:height(traveling_ERP_table)
            direction = traveling_ERP_table{row, "Direction"};
            alignment = traveling_ERP_table{row, "Alignment"};
            contacts = traveling_ERP_table{row, "Contacts"};
            ct_st = cur_elec_contact_ind(contacts(1)); 
            ct_end = cur_elec_contact_ind(contacts(2));
            subwin_st = traveling_ERP_table{row, 3};
            subwin_end = traveling_ERP_table{row,4};

            data_out_dir = sprintf('%s/data_filtered_trial_avg/%s',out_dir, alignment);
            % Record max amplitude (of all contacts participating in the wave)
            max_amp = 0;
            for wave_ct = ct_st:ct_end
                % load real-valued trial avg data for wave_ct 
                name = elec_name{wave_ct};
                fname_data_avg = sprintf("%s/realavg_%s_%s_ch%03d.mat", data_out_dir, alignment, name, wave_ct);
                load(fname_data_avg, "avg_events")
                cur_max = max(abs(avg_events(subwin_st:subwin_end)));
                if cur_max > max_amp
                    max_amp = cur_max;
                end
            end
            max_amps = [max_amps; max_amp]; % add a column to the trav wave table for max_amp
           
        end % end loop through rows of trav ERP table
        traveling_ERP_table.MaxAmp = max_amps;
        fname_table = sprintf('%s/%s_trialavg_travERP_table.mat', table_out_dir, cur_letter); 
        save(fname_table, "traveling_ERP_table")

    end % end loop through electrodes

end % end loop through subj