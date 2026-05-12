%% Plot single trial data - first four unique inspections

subject_ID = 'EMU038';
reference = 'Ground'; Fs = 2048; num_sessions = 1; sesnum = 1;
win = [0, Fs/2];

alignments = ["first_unique_attribute","second_unique_attribute","third_unique_attribute","fourth_unique_attribute"];
align_names = ["first attribute","second attribute","third attribute","fourth attribute"];

load(sprintf('%s/%s/%s_MAD_SES%d_Setup.mat',data_base_dir,subject_ID,subject_ID,1));

% NOTE: assuming only one recording session         
% [align_times,trial_numbers] = get_align_times(filters, trial_times, trial_words, 'inspection');
% remove_ind = isnan(align_times);
% align_times(isnan(align_times)) = [];
% trial_numbers(remove_ind) = [];
% align_times = round(align_times*Fs); 
% num_trials = max(trial_numbers); % number of trials this patient played

for align_num = 1:4
    [align_times,trial_numbers] = get_align_times(filters, trial_times, trial_words, alignments(align_num));
    remove_ind = isnan(align_times);
    align_times(isnan(align_times)) = [];
    trial_numbers(remove_ind) = [];
    align_times = round(align_times*Fs);
    eval(sprintf('align_times%d = [align_times trial_numbers];',align_num))

    % if strcmp(subject_ID,'EMU001')
    %     [align_times_inspect,trial_numbers_inspect] = get_align_times(filters, trial_times, trial_words, '2opt_2att_inspection');
    % 
    %     remove_ind = isnan(align_times_inspect);
    %     align_times_inspect(isnan(align_times_inspect)) = [];
    %     trial_numbers_inspect(remove_ind) = [];
    % 
    %     align_times_inspect = round(align_times_inspect*Fs); 
    % 
    %     align_times = intersect(align_times, align_times_inspect);
    % end
end
all_trials = unique(align_times4(:,2));
if numel(all_trials) ~= numel(align_times4(:,2))
    disp("Warning: duplicate trials for fourth unique attribute inspection")
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


for e = 7 %1:numel(elec_letters)  TODO: [6 7 8]
    cur_letter = elec_letters(e);
    
    for trial = all_trials(103:end)' %all_trials'

        fg1 = figure; 
        set( fg1, 'position', [ 88  100  1450  760 ] ) 
        tiledlayout(2,2,'TileSpacing','compact');
        %title(t,sprintf('Trial %d',trial))
    
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
        colors = [
                    0.76 0.86 0.98;
                    0.72 0.83 0.96;
                    0.68 0.79 0.94;
                    0.64 0.76 0.92;
                    0.60 0.72 0.90;
                    0.56 0.69 0.88;
                    0.52 0.65 0.86;
                    0.48 0.62 0.84;
                    0.44 0.58 0.82;
                    0.40 0.55 0.80;
                    0.36 0.51 0.78;
                    0.32 0.48 0.76;
                    0.28 0.44 0.74;
                    0.24 0.41 0.72;
                    0.20 0.37 0.70;
                    0.16 0.34 0.68;
                    0.12 0.30 0.66;
                    0.06 0.26 0.63];
    
        for inspect_num = 1:4
            eval(sprintf('align_times_trials = align_times%d;',inspect_num))
            cur_trial_idx = find(align_times_trials(:,2) == trial);
            if isempty(cur_trial_idx), close; break; end
            if numel(cur_trial_idx) > 1, fprintf("WARNING: more than one unique inspection %d associated with trial %d", inspect_num, trial); end
            align_time = align_times_trials(cur_trial_idx,1);
    
            out_dir = sprintf('/media/Data/Human_Intracranial_MAD/analysis/Cross_Freq_Coupling/%s/%s/synchrony_data/generalized_phase',subject_ID, reference);  
            final_out_dir = sprintf('%s/waveform_plots/single_trial/%s',out_dir, 'inspection');
            %final_out_dir = sprintf('%s/waveform_plots/trial_averaged/%s',out_dir, alignment);
            if ~exist(final_out_dir, 'dir'), mkdir(final_out_dir); end
        
            data_out_dir = sprintf('%s/data_filtered_trial_avg/%s',out_dir, alignment);
            if ~exist(data_out_dir, 'dir'), mkdir(data_out_dir); end
    
            nexttile;
            for cnum = 1:numel(contacts) 
                contact = contacts(cnum);
                cur_name = cur_elec_contact_names(cnum);
                
                    
                load(sprintf('%s/%s/separate_channel_files/%s_MAD_SES%d_ch%03d.mat',data_base_dir,subject_ID,subject_ID,sesnum,contact)); % for saving phase for ALL contacts

                % broadband filter
                [b,a] = butter( 4, [5 40] ./ (Fs/2) ); filtered = filtfilt( b, a, data ); 
                data = filtered;
                data_event = data((round(align_time-win(1)):round(align_time+win(2)-1)));
                data_event = ( data_event - mean(data_event) ) ./ std(data_event);     

                % Construct complex signal
                lp = 5;
                xgp_avg = generalized_phase_vector( avg_events, Fs, lp );
        
                % PLOT AVERAGED TIME SERIES
                plot((1:numel(data_event))/Fs,data_event,'linewidth',2,'color',colors(cnum,:)) % -1024:-1
                hold on
    
                xlim([0,0.5]); 
                %xlim([0,0.3]); % TEMP 12/19/25
                ylim([-4, 4]);
                xlabel('time from alignment (s)'); 
                ylabel('zscored amplitude');
                %yl1 = ylim;
                title(sprintf("Trial %d %s",trial, align_names(inspect_num)))
                %title(sprintf("%s Elec %s TimeSeries",align_names(align_num),cur_letter))
                
            end % end loop through contacts
        end % end loop through alignments
        fname = sprintf("%s/elec%s_trial%d_inspections.jpg",final_out_dir, cur_letter, trial);
    end % end loop through trials
end