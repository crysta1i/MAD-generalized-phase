%% Induced waves

%% 1. Trial-average after aligning all single-trial signals to highest peak
% Find the time of the highest peak/lowest trough (of all contacts on a probe) 
% in the trial-averaged signal. (The contact that reaches that highest peak/trough 
% is the reference contact.) For each single-trial signal on the reference contact, 
% find the peak/trough closest to the highest one (as identified in the trial-averaged 
% signal) and align all single-trial signals by that peak/trough. Record the time 
% offsets that produce these re-alignments and reuse them on all other contacts on 
% the probe.
subject_ID = 'EMU038'; reference = 'Ground';

% Pre-requisites: 
% - save the trial-average (real-valued) signals for a given alignment
% - save the full time series xgp signals for all contacts (extract single events from this)
datapath = '/media/Data/Human_Intracranial_MAD/';
data_base_dir = sprintf('%s1_formatted',datapath); 
out_dir = sprintf('/media/Data/Human_Intracranial_MAD/analysis/Cross_Freq_Coupling/%s/%s/synchrony_data/generalized_phase',subject_ID, reference);          

alignments = ["inspection","single_opt_first_inspection","full_single_opt_info"];
align_names = ["Inspection","Single Opt 1st Inspection","Full Single Opt Info"];
for align_num = 1:numel(alignments)
    alignment = alignments(align_num);
    final_out_dir = sprintf('%s/waveform_plots/trial_averaged/%s',out_dir, alignment);
    if ~exist(final_out_dir, 'dir'), mkdir(final_out_dir); end
    
    data_out_dir = sprintf('%s/data_filtered_trial_avg/%s',out_dir, alignment);
    load(sprintf('%s/%s/%s_MAD_SES%d_Setup.mat',data_base_dir,subject_ID,subject_ID,1));
    
    % get elec letters
    elec_letters = [];
    channel_names_bs = elec_name;
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


        max_amp_all = 0; max_amp_time = 0; max_amp_ct = '';
        min_amp_all = 0; min_amp_time = 0; min_amp_ct = '';
        for cnum = 1:numel(cur_elec_contact_names)
            cur_name = cur_elec_contact_names(cnum);
            contact = cur_elec_contact_ind(cnum);
            
            fname_data_avg = sprintf("%s/realavg_%s_%s_ch%03d.mat", data_out_dir, alignment, cur_name, contact);
            load(fname_data_avg)
        
            [max_amp, max_time] = max(avg_events); % for now, just try aligning to highest peak and lowest trough
            [min_amp, min_time] = min(avg_events);
            if max_amp > max_amp_all
                max_amp_all = max_amp; max_amp_time = max_time; 
                max_amp_ct = cur_elec_contact_names(cnum); max_amp_ct_ind = cur_elec_contact_ind(cnum);
            end
            if min_amp < min_amp_all
                min_amp_all = min_amp; min_amp_time = min_time; 
                min_amp_ct = cur_elec_contact_names(cnum); min_amp_ct_ind = cur_elec_contact_ind(cnum);
            end
        end % end loop through contacts on the electrode
        
        % NOTE: min_amp_ct and max_amp_ct will likely be different!!
        
        
        peak_align_times = []; %trough_align_times = [];
        peak_align_data = []; %trough_align_data = [];
        for sesnum = 1:num_sessions
            load(sprintf('%s/%s/%s_MAD_SES%d_Setup.mat',data_base_dir,subject_ID,subject_ID,sesnum));
            
            [align_times,trial_numbers] = get_align_times(filters, trial_times, trial_words, alignment);
    
            remove_ind = isnan(align_times);
            align_times(isnan(align_times)) = [];
            trial_numbers(remove_ind) = [];
            
            align_times = round(align_times*Fs);  
    
            cur_name = max_amp_ct; contact = max_amp_ct_ind;
            
            % load real-valued signal for the REFERENCE contact
            load(sprintf('%s/%s/separate_channel_files/%s_MAD_SES%d_ch%03d.mat',data_base_dir,subject_ID,subject_ID,sesnum,contact)); % for saving phase for ALL contacts
            [b,a] = butter( 4, [5 40] ./ (Fs/2) ); data = filtfilt( b, a, data );

            % load xgp for the REFERENCE contact
            load(sprintf('%s/extracted_phase/xgp_%s_ses%02d_ch%03d.mat',out_dir,cur_name,sesnum,contact));
            
            closest_pk_times = zeros(numel(align_times),1); 
            %closest_tgh_times = zeros(numel(align_times),1);
            closest_pk_data = zeros(Fs/2,numel(align_times));
            %closest_tgh_data = zeros(Fs/2,numel(align_times));
            for event = 1:numel(align_times)  
                data_event = data((round(align_times(event)-win(1)):round(align_times(event)+win(2)-1)));
                data_event = ( data_event - mean(data_event) ) ./ std(data_event);

                xgp_event = xgp((round(align_times(event)-win(1)):round(align_times(event)+win(2)-1)));
                % in the single trial data, find all peaks/troughs, then identify the time of the peak closest to the one identified as max_amp
                posneg_xgp = angle(clipped_xgp);
                posneg_xgp(posneg_xgp >= 0) = 1; posneg_xgp(posneg_xgp < 0) = -1;
                peaks_xgp = find(diff(posneg_xgp) > 0); %troughs_xgp = find(diff(posneg_xgp) < 0);
                [~,idx_min_peak] = min(abs(peaks_xgp - max_time)); closest_pk_time = peaks_xgp(idx_min_peak);
                %[~,idx_min_trough] = min(abs(peaks_xgp - min_time)); closest_tgh_time = troughs_xgp(idx_min_trough);
                
                % need to save closest_pk_times/closest_tgh_times to use on other contacts
                closest_pk_times(event) = closest_pk_time; %closest_tgh_times(event) = closest_tgh_time;

                start_time_peaks = round(align_times(ii)-win(1)) + closest_pk_time - 1 - 100;
                %start_time_troughs = round(align_times(ii)-win(1)) + closest_tgh_time - 1 - 100;
                cur_event_data_pk = data(start_time_peaks:start_time_peaks + abs(win(1)-win(2)) - 1);
                %cur_event_data_tgh = data(start_time_troughs:start_time_troughs + abs(win(1)-win(2)) - 1);
                closest_pk_data(:,event) = (cur_event_data_pk - mean(cur_event_data_pk)) ./ std(cur_event_data_pk);
                %closest_tgh_data(:,event) = (cur_event_data_tgh - mean(cur_event_data_tgh)) ./ std(cur_event_data_tgh);
                    % start 100 samples BEFORE alignment to peak/trough?
                    % end the event window Fs/2 samples after alignment

                
            end % end loop through events
            peak_align_data = [peak_align_data closest_pk_data];
            %trough_align_data = [trough_align_data closest_tgh_data];
            peak_align_times = [peak_align_times; closest_pk_times];
            %trough_align_times = [trough_align_times; closest_tgh_times];
        end % end loop through sessions
        % average peak_align_data and trough_align_data horizontally
        peak_realigned_avg = mean(peak_align_data,2);
        %trough_realigned_avg = mean(trough_align_data,2);

        all_ct_peak_realigned = zeros(abs(win(1)-win(2)),numel(cur_elec_contact_ind));
        for cnum = 1:numel(cur_elec_contact_ind)
            cur_ct_peak_realigned = [];
            tot_num_events = 0;
            contact = cur_elec_contact_ind(cnum);
            for sesnum = 1:num_sessions
                
                load(sprintf('%s/%s/%s_MAD_SES%d_Setup.mat',data_base_dir,subject_ID,subject_ID,sesnum));
            
                [align_times,trial_numbers] = get_align_times(filters, trial_times, trial_words, alignment);
        
                remove_ind = isnan(align_times);
                align_times(isnan(align_times)) = [];
                trial_numbers(remove_ind) = [];
                
                align_times = round(align_times*Fs);

                load(sprintf('%s/%s/separate_channel_files/%s_MAD_SES%d_ch%03d.mat',data_base_dir,subject_ID,subject_ID,sesnum,contact)); 
                [b,a] = butter( 4, [5 40] ./ (Fs/2) ); data = filtfilt( b, a, data );
                
                cur_ct_ses_peak_realigned = zeros(abs(win(1)-win(2)),numel(align_times));
                for event = 1:numel(align_times)
                    
                    cur_event_start = peak_align_times(tot_num_events + event);
                    data_event = data(cur_event_start:cur_event_start + abs(win(1)-win(2)) - 1);
                    data_event = ( data_event - mean(data_event) ) ./ std(data_event);
            
                    % tot_num_events+1:tot_num_events+numel(align_times)
                    cur_ct_ses_peak_realigned(:,event) = data_event;
                end
                cur_ct_peak_realigned = [cur_ct_peak_realigned cur_ct_ses_peak_realigned];
                tot_num_events = tot_num_events + numel(align_times);
            end % end loop through sessions
            cur_ct_peak_realigned = mean(cur_ct_peak_realigned,2);
            all_ct_peak_realigned(:,cnum) = cur_ct_peak_realigned;

            
        end % end loop through contacts

        % PLOT
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
        fg1 = figure; 
        set( fg1, 'position', [ 88  100  1400  600 ] ) %[ 88  1593  1250  420 ]
        for cnum = 1:numel(cur_elec_contact_ind)
            plot(all_ct_peak_realigned(:,cnum),'linewidth',2,'color',colors(cnum,:)) % -1024:-1
            hold on
            xlim([0,abs(win(1)-win(2))]); 
            %ylim([-0.5, 0.5]);
            xlabel('samples from alignment'); 
            ylabel('zscored amplitude');
            yl1 = ylim;
            title(sprintf("%s Elec %s Peak-Aligned TimeSeries",align_names(align_num),cur_letter))
        end

        fname = sprintf("%s/align2peak_elec%s.jpg", final_out_dir,cur_letter);
        print(fg1,'-djpeg',fname)
        fprintf("saved %s elec %s \n",alignment,cur_letter)
        close;

    end % end loop through electrodes
end % end loop through alignments