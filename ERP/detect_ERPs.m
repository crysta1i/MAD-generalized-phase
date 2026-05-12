%% Statistical detection of ERPs + FDR correction

function ERP_table = detect_ERPs(subject_ID, reference, alignments, cur_letter)

% 3/4/26: This is the version of the script in TravWaves/Code/ERP
% Latest edited version

% Inputs
% - subject_ID, reference, alignment
% - probe letter (e.g., "Z'", "F", "N")

% Outputs
% to save: ERP_table (one per contact and alignment)
% - Columns: peak/trough, #samples (time) where significant ERP occurs, idx-before, idx-after, amplitude, p-value

% Permutation testing (Compare aligned to time-jittered time series)
addpath('/media/Data/Human_Intracranial_MAD/analysis/TravWaves/Code')

[~, data_base_dir, tw_out_dir, ~, num_sessions, Fs, ~] = tw_setup(subject_ID, reference);
win = [0 Fs/2];

%load(sprintf('%s/%s/%s_MAD_SES%d_Setup.mat',data_base_dir,subject_ID,subject_ID,1));

[cur_elec_contact_ind, cur_elec_contact_names] = get_single_probe_contacts(reference, subject_ID, cur_letter);

for alignment = alignments

    for cnum = 1:numel(cur_elec_contact_ind) 
        contact = cur_elec_contact_ind(cnum);
        cur_name = cur_elec_contact_names(cnum);
        
        num_events = 0;
        sum_events_padded = zeros(max(win) + 2*Fs,1);
        for sesnum = 1:num_sessions
            load(sprintf('%s/%s/%s_MAD_SES%d_Setup.mat',data_base_dir,subject_ID,subject_ID,sesnum), "filters", "trial_times", "trial_words");
            
            [align_times,~] = get_align_times(filters, trial_times, trial_words, alignment);
            align_times(isnan(align_times)) = [];
            align_times = round(align_times*Fs);  

            if strcmp(subject_ID,'EMU001')
                [align_times_inspect,~] = get_align_times(filters, trial_times, trial_words, '2opt_2att_inspection');
                align_times_inspect(isnan(align_times_inspect)) = [];

                align_times_inspect = round(align_times_inspect*Fs); 

                align_times = intersect(align_times, align_times_inspect);
            end
            
            load(sprintf('%s/%s/separate_channel_files/%s_MAD_SES%d_ch%03d.mat',data_base_dir,subject_ID,subject_ID,sesnum,contact), "data");
            
            % broadband filter
            [b,a] = butter( 4, [5 40] ./ (Fs/2) ); filtered = filtfilt( b, a, data ); 
            data = filtered; 

            for ii = 1:numel(align_times)

                %data_event = data((round(align_times(ii)-win(1)):round(align_times(ii)+win(2)-1)));
                data_event_padded = data(align_times(ii)-win(1)-Fs:align_times(ii)+win(2)-1+Fs);
                data_event_padded = ( data_event_padded - mean(data_event_padded) ) ./ std(data_event_padded);
                sum_events_padded = sum_events_padded + data_event_padded;
    
            end % end loop through events
            num_events = num_events + numel(align_times);
        end % end loop through sessions

        avg_events_aligned_padded = sum_events_padded ./ num_events;
        avg_events_aligned = avg_events_aligned_padded(1+Fs:1+Fs+win(2)-1);

        % FIND ALL PEAKS AND TROUGHS
        % Construct complex signal
        lp = 5; xgp_avg = generalized_phase_vector( avg_events_aligned_padded, Fs, lp );  
        xgp_avg = xgp_avg(1+Fs:1+Fs+win(2)-1);
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
        peak_trough_distr = time_jitter(subject_ID, reference, iterations, contact, alignment);
        peak_trough_distr = sort(peak_trough_distr,'ascend');
        
        pv_table = cell(0,6); 
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
            
            pv_table(end+1,:) = {"peak", peak_idx, idx_trough_before, idx_trough_after, peak, pv};
        end
            
        for trough_num = 1:numel(troughs_align)
            trough = troughs_align(trough_num);
            %if trough < peak_trough_distr(critical_trough_idx)
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
                %ERP_table(end+1,:) = {alignment, cur_letter, cur_name, contact, "trough", trough_idx, trough, idx_peak_before, idx_peak_after,pv};
                pv_table(end+1,:) = {"trough", trough_idx, idx_peak_before, idx_peak_after, trough, pv};
            %end
        end
        pv_table = cell2table(pv_table, "VariableNames",["Peak/Trough","Idx","Idx Before","Idx After","Amplitude","p-value"]);

        % FDR CORRECTION
        pvals = pv_table(:,"p-value");
        pvals = table2array(pvals);

        q = 0.05; method = 'dep'; report = 'no';
        [h, ~, ~, adj_p]=fdr_bh(pvals,q,method,report);

        ERP_table = addvars(pv_table, adj_p, 'NewVariableNames', 'Adj p-val');

        ERP_table(h == 0, :) = []; % remove rows of ERP table where h == 0 (not sig)

        if ~exist(sprintf('%s/ERP/%s',tw_out_dir,alignment), 'dir'), mkdir(sprintf('%s/ERP/%s',tw_out_dir,alignment)); end
        fname_table = sprintf('%s/ERP/%s/sig_ERPs_table_%s_%03d.mat', tw_out_dir,alignment,cur_name,contact); 
        save(fname_table, "ERP_table")
        %disp("table saved")

    end % end loop through contacts

end % end loop through alignments

%% TEMP: test xgp vs derivative method for estimating peak/trough times
% % CONCLUSION from this: use angle(xgp) to estimate peak times, not real
% % valued signal
% 
% contact = 3; subject_ID = 'EMU038'; cur_name = "P'left13"; reference = 'Ground';
% 
% % xgp
% out_dir = sprintf('/media/Data/Human_Intracranial_MAD/analysis/Cross_Freq_Coupling/%s/%s/synchrony_data/generalized_phase',subject_ID, reference);
% fname = sprintf('%s/extracted_phase/xgp_%s_ses%02d_ch%03d.mat',out_dir,cur_name,1,contact);
% load(fname) % loads xgp
% clipped_xgp = xgp(1:2000); posneg_xgp = angle(clipped_xgp);
% posneg_xgp(posneg_xgp >= 0) = 1; posneg_xgp(posneg_xgp < 0) = -1;
% peaks_xgp = find(diff(posneg_xgp) > 0); troughs_xgp = find(diff(posneg_xgp) < 0);
%     % negative to positive phase (cross 0) = peak
%     % positive to negative phase (+pi to -pi) = trough
% 
% % derivative
% datapath = '/media/Data/Human_Intracranial_MAD/'; data_base_dir = sprintf('%s1_formatted',datapath);
% load(sprintf('%s/%s/separate_channel_files/%s_MAD_SES%d_ch%03d.mat',data_base_dir,subject_ID,subject_ID,1,contact)); % for saving phase for ALL contacts
% % filter the data from 5-40 Hz
% [b,a] = butter( 4, [5 40] ./ (Fs/2) ); data = filtfilt( b, a, data );
% clipped_real = data(1:2000); deriv_real = diff(clipped_real);
% posneg_deriv = deriv_real; 
% posneg_deriv(posneg_deriv > 0) = 1; % increasing
% posneg_deriv(posneg_deriv < 0) = -1; % decreasing
% peaks_real = find(diff(posneg_deriv) < 0) + 1; troughs_real = find(diff(posneg_deriv > 0)) + 1; 

end