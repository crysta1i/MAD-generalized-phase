% Compare ground and neighbor-average referenced single-trial data - filtered (5-40Hz)
% Run this after running detect_waves_trialavg.m -- 1.25.26 section
    % need to set up probe letter, contacts_st/end, subwin st/end

references = {'Ground','neighbor_average'}; ref_names = {'Ground','neighbor average'};
alignment = "fourth_unique_attribute"; align_name = 'fourth unique attribute';
cur_letter = 'N';
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
    
            % 2/21/26 -- no longer using angle(xgp) method to find peaks/troughs
    
            % xgp = load(sprintf('%s/extracted_phase/xgp_%s_ses%02d_ch%03d.mat',out_dir,cur_name,sesnum,contact));
            % xgp = xgp.xgp;
            % xgp_event = xgp(align_times(rand_event)-win(1):align_times(rand_event)+win(2)-1);
            % xgp_event_subwin = xgp_event(subwin_st:subwin_end);
            % angle_subwin = angle(xgp_event_subwin);
    
            cur_ct_pk_times = [];
            if peak_wave
                cur_ct_pk_times = find(sign(diff([data_event; data_event(end)]))<0&sign(diff([data_event(1); data_event]))>0);
            else
                cur_ct_pk_times = find(sign(diff([data_event; data_event(end)]))>0&sign(diff([data_event(1); data_event]))<0);
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