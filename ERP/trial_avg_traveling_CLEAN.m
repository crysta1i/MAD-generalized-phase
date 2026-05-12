% trial_avg_travling_CLEAN.m

% Options

%reference = 'Ground';
reference = 'neighbor_average';

alignments = ["full_single_opt_info"]; align_names = ["Full Single Option Info"];
%alignments = ["inspection"]; align_names = ["Inspection"];
%alignments = ["probability","amount"]; align_names = ["Probability","Amount"];
%alignments = ["likely_outcome","unlikely_outcome"];
%alignments = ["positive", "negative"];
% alignments = ["single_opt_first_inspection","full_single_opt_info"]; align_names = ["Single Option First Inspection","Full Single Option Info"];
% alignments = ["first_unique_attribute", "second_unique_attribute","third_unique_attribute","fourth_unique_attribute"]; align_names = ["Option 1 First Inspection", "Option 1 Second Inspection","Option 2 First Inspection", "Option 2 Second Inspection"];

colors = [  0.76 0.86 0.98;
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

addpath('/media/Data/Human_Intracranial_MAD/_toolbox/spectral-analysis');
addpath('/media/Data/Human_Intracranial_MAD/analysis/PAC_code/CF_Coupling/generalized phase/generalized-phase/analysis/')
subject_IDs = {'EMU001','EMU024', 'EMU025','EMU030','EMU037','EMU038','EMU039','EMU041','EMU047','EMU051'}; 

datapath = '/media/Data/Human_Intracranial_MAD/';
data_base_dir = sprintf('%s1_formatted',datapath); 

for subject_num = 6 %[1 4 5 7 8 9 10] %6 %4:numel(subject_IDs)
    
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

    if strcmp(reference,'Ground')
        load(sprintf('%s/%s/%s_MAD_SES%d_Setup.mat',data_base_dir,subject_ID,subject_ID,1));
        channel_ind = (1:numel(elec_name))';
    else
        load(sprintf('%s/%s/%s_MAD_SES%d_Setup_%s.mat',data_base_dir,subject_ID,subject_ID,1,reference));
        elec_name = channel_name; elec_area = channel_area; elec_region = channel_region; elec_location = channel_location;
    end

    % get elec letters
    elec_letters = [];
    for i = 1:numel(elec_name) 
        cur_name = elec_name{i};
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

    for e = 1 %:numel(elec_letters) 
        cur_letter = 'N';% TEMP 
        %cur_letter = elec_letters(e);

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
            cur_elec_contact_ind = [cur_elec_contact_ind; channel_ind(j)];
            
            j = j+1; 
            if j > numel(elec_name), break; end
            cur_name = elec_name{j};
        end

        % Ground: lateral contacts listed first (e.g., P15 -> P1)
        % Neighbor-average: medial contacts listed first (P1 -> P15)
        if strcmp(reference,'neighbor_average')
            cur_elec_contact_names = flip(cur_elec_contact_names);
            cur_elec_contact_ind = flip(cur_elec_contact_ind);
        end
 
        for align_num = 1:numel(alignments)
            alignment = alignments(align_num);

            out_dir = sprintf('/media/Data/Human_Intracranial_MAD/analysis/Cross_Freq_Coupling/%s/%s/synchrony_data/generalized_phase',subject_ID, reference);  
            plot_out_dir = sprintf('%s/waveform_plots/trial_averaged/%s',out_dir, 'inspection');
            if ~exist(plot_out_dir, 'dir'), mkdir(plot_out_dir); end
        
            data_out_dir = sprintf('%s/data_filtered_trial_avg/%s',out_dir, alignment);
            if ~exist(data_out_dir, 'dir'), mkdir(data_out_dir); end
            
            to_plot_aligned = [];
            to_plot_jittered = [];
            for cnum = 1:numel(cur_elec_contact_ind) 
                contact = cur_elec_contact_ind(cnum);
                cur_name = cur_elec_contact_names(cnum);
                
                sum_events = zeros(max(win),1);
                sum_events_off = zeros(max(win),1);
                num_events = 0;
                for sesnum = 1:num_sessions      
                    if strcmp(reference,'Ground')
                        load(sprintf('%s/%s/%s_MAD_SES%d_Setup.mat',data_base_dir,subject_ID,subject_ID,sesnum));
                    else
                        load(sprintf('%s/%s/%s_MAD_SES%d_Setup_%s.mat',data_base_dir,subject_ID,subject_ID,sesnum,reference));
                        elec_name = channel_name; elec_area = channel_area; elec_region = channel_region; elec_location = channel_location;
                    end

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
                    
                    if strcmp(reference,'Ground')
                        load(sprintf('%s/%s/separate_channel_files/%s_MAD_SES%d_ch%03d.mat',data_base_dir,subject_ID,subject_ID,sesnum,contact)); 
                    elseif strcmp(reference,'neighbor_average')
                        load(sprintf('%s/%s/separate_channel_files/%s/%s_MAD_SES%d_ch%03d.mat',data_base_dir,subject_ID,reference,subject_ID,sesnum,contact));
                    end
                    
            
                    % % SUBTRACT MOVING AVERAGE
                    % [n_ts, n_channels] = size(data); % n_ts = number of time steps total
                    % win_len = round(mean_subtract_window*Fs);
                    % ma = zeros(size(data));
                    % for i = 1:n_channels
                    %     ma(:,i) = movmean(data(:,i), win_len);data(:,i);
                    % end
                    % mean_subtracted = data - ma;
        
                    % %FILTER 60 HZ NOISE
                    % d = designfilt('bandstopiir','FilterOrder',2, ...
                    %     'HalfPowerFrequency1',59,'HalfPowerFrequency2',61, ...
                    %     'DesignMethod','butter','SampleRate',Fs);
                    % d1 = designfilt('bandstopiir','FilterOrder',2, ...
                    %     'HalfPowerFrequency1',119,'HalfPowerFrequency2',121, ...
                    %     'DesignMethod','butter','SampleRate',Fs);
                    % %filtered = filtfilt(d, mean_subtracted);
                    % filtered = filtfilt(d, data); % TEMP
                    % filtered = filtfilt(d1, filtered);
                    % data = filtered; % TEMP
                    
                    % broadband filter
                    [b,a] = butter( 4, [5 40] ./ (Fs/2) ); filtered = filtfilt( b, a, data ); 
                    data = filtered;
  
                    num_events = num_events + numel(align_times);
                    for ii = 1:numel(align_times) 
        
                        data_event = data((round(align_times(ii)-win(1)):round(align_times(ii)+win(2)-1)));
                        data_event = ( data_event - mean(data_event) ) ./ std(data_event);
                        sum_events = sum_events + data_event;      
            
                    end

                    % RNG for random time offsets
                    rng(10, 'twister') % seed = 10
                    rand_nums = rand(numel(align_times),1);

                    for ii = 1:numel(align_times)

                        offset_sec = -0.25 + 0.5 * rand_nums(ii);
                        offset_samples = offset_sec * Fs;

                        data_event_off = data((round(align_times(ii)-win(1)+offset_samples):round(align_times(ii)+win(2)-1+offset_samples)));
                        data_event_off = (data_event_off - mean(data_event_off))./ std(data_event_off); % zscore - TEMP
                        sum_events_off = sum_events_off + data_event_off;     

                    end
                    

                end % end loop through sessions

                avg_events = sum_events./ num_events;
                avg_events_off = sum_events_off./ num_events;
                to_plot_aligned = [to_plot_aligned avg_events];
                to_plot_jittered = [to_plot_jittered avg_events_off];

                % Construct complex signal
                lp = 5;
                xgp_avg = generalized_phase_vector( avg_events, Fs, lp );
    
                fname_xgp_avg = sprintf("%s/centered_%s_xgp_%s_ch%03d.mat", data_out_dir, alignment, cur_name, contact);
                save(fname_xgp_avg, "xgp_avg")
                
                fname_data_avg = sprintf("%s/realavg_%s_%s_ch%03d.mat", data_out_dir, alignment, cur_name, contact);
                save(fname_data_avg, "avg_events")
                
            end % end loop through contacts

            % PLOT
            fg1 = figure; 
            set( fg1, 'position', [ 88  100  1450  760 ] )
            tiledlayout(2,1,'TileSpacing','compact');

            nexttile;
            for cnum = 1:width(to_plot_aligned)
                plot((1:height(to_plot_aligned))/Fs,to_plot_aligned(:,cnum),'linewidth',2,'color',colors(cnum,:)) % -1024:-1
                hold on
    
                xlim([0,height(to_plot_aligned)/Fs]);
                %xlim([0,Fs/2]); 
                %ylim([-0.5, 0.6]);
                xlabel('time from alignment (s)'); 
                ylabel('zscored amplitude');
                yl1 = ylim;
                title(sprintf("%s",align_names(align_num)))
            end

            nexttile;
            for cnum = 1:width(to_plot_jittered)
                plot((1:height(to_plot_jittered))/Fs,to_plot_jittered(:,cnum),'linewidth',2,'color',colors(cnum,:)) % -1024:-1
                hold on
    
                xlim([0,height(to_plot_jittered)/Fs]); 
                %ylim([-0.5, 0.6]);
                xlabel('time from alignment (s)'); 
                ylabel('zscored amplitude');
                ylim(yl1);
                title(sprintf("%s time-jittered",align_names(align_num)))
            end
        
        end % end loop through alignments
        %fname = sprintf("%s/allses_%s_vs_%s_elec%s.jpg", final_out_dir, alignments(1), alignments(2),cur_letter);
        fname = sprintf("%s/allses_%s_elec%s.jpg", plot_out_dir, alignments(1),cur_letter);
        %fname = sprintf("%s/all4inspections_elec%s.jpg", final_out_dir,cur_letter);
        print(fg1,'-djpeg',fname)
        close;
    end % end loop through electrodes

end %end loops through subj