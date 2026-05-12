% Intersite Phase Clustering / PLV analysis. Updated 3/17/26


subject_ID = 'EMU038'; reference = 'neighbor_average';
alignments = ["first_unique_attribute", "second_unique_attribute","third_unique_attribute","fourth_unique_attribute"]; 

cur_letter = 'N';
[~, ~, tw_out_dir, ~, ~, Fs, ~] = tw_setup(subject_ID, reference);
plot_out_dir = sprintf("%s/Plots/%s/%s/ISPC", tw_out_dir, reference, subject_ID);
if ~exist(plot_out_dir, 'dir'), mkdir(plot_out_dir); end
[cur_elec_contact_ind, cur_elec_contact_names] = get_single_probe_contacts(reference, subject_ID, cur_letter);


% ------------------------- 1. TRIAL-AVG ISPC -------------------------------

% for alignment = alignments
%     align_name = strrep(alignment, '_', ' '); 
    
%     adj_pairs_ISPC_avg = []; % each row = time series of PLV for one pair of contacts
%     adj_pairs_lag_avg = [];
%     pair_labels = strings(1, numel(cur_elec_contact_ind)-1);
%     for cnum = 1:numel(cur_elec_contact_ind)-1 % for each adjacent pair of contacts
%         [ISPC_avg, lag_avg] = ISPC_across_trials(subject_ID, reference, alignment, cur_letter, cnum, cnum+1);
%         adj_pairs_ISPC_avg = [adj_pairs_ISPC_avg; ISPC_avg.'];
%         adj_pairs_lag_avg = [adj_pairs_lag_avg; lag_avg.'];
%         pair_labels(cnum) = sprintf("%s-%s", cur_elec_contact_names(cnum), cur_elec_contact_names(cnum+1));
%     end
%     time = (1:width(adj_pairs_lag_avg))/Fs;

%     % --------- Heat maps ------------
%     % figure;
%     % imagesc(time, 1:height(adj_pairs_ISPC_avg), adj_pairs_ISPC_avg); 
%     % colormap hot                   
%     % c = colorbar;
%     % c.Label.String = 'PLV'; 
%     % clim([0 1]);
%     % xlabel('Time (s)');
%     % ylabel('Electrode pair');
%     % title(sprintf("PLV/ISPC adjacent contacts probe %s, %s", cur_letter, align_name));

%     % yticks(1:numel(cur_elec_contact_ind)-1);
%     % set(gca, 'YTickLabel', pair_labels(yticks)); 

%     % filename = sprintf("%s/probe%s_%s_avg_heatmap.jpg", plot_out_dir, cur_letter, alignment);
%     % print(gcf, filename, '-djpeg')
%     % close;

%     % --------- Line plots ------------
%     M = size(adj_pairs_lag_avg,1); N = size(adj_pairs_lag_avg,2);
%     % fg1 = figure('Units','normalized','Position',[0.1 0.1 0.6 0.8]);
%     % ax = gobjects(M,1);
%     % top = 0.95; bottom = 0.06; left = 0.12; right = 0.92;
%     % totalHeight = top - bottom;
%     % h = totalHeight / M;          % height per y-axis

%     % for k = 1:M
%     %     % compute vertical position; want row 1 at top
%     %     pos = [left, top - k*h, right-left, h*0.95]; % spacing
%     %     ax(k) = axes('Position', pos);
%     %     plot(time, adj_pairs_ISPC_avg(k,:), 'k');
%     %     xlim([time(1) time(end)]);
%     %     ylim([min(adj_pairs_ISPC_avg(k,:)) max(adj_pairs_ISPC_avg(k,:))]); % use per-row limits
%     %     set(ax(k), 'YTick', []); 
%     %     % label each row on left
%     %     text(time(1) - 0.02*(time(end)-time(1)), mean(ylim), pair_labels(k), 'HorizontalAlignment','right', 'VerticalAlignment','middle');
%     %     box off;
%     % end
%     % linkaxes(ax, 'x');
%     % xlabel('Time (s)');
%     % fname = sprintf("%s/probe%s_%s_avg_lineplots.jpg", plot_out_dir, cur_letter, alignment);
%     % print(fg1,'-djpeg',fname)
%     % close;


%     % --------- ISPC + lag subplots ------------
%     fg1 = figure('Units','normalized','Position',[0.1 0.1 0.8 0.5]); % subplots (left=ISPC as heatmap, right=lag as lineplot)
%     tiledlayout(1,2,'TileSpacing','compact');
    
%     nexttile;
%     imagesc(time, 1:height(adj_pairs_ISPC_avg), adj_pairs_ISPC_avg); 
%     colormap hot                   
%     c = colorbar;
%     c.Label.String = 'PLV'; 
%     clim([0 1]);
%     xlabel('Time (s)');
%     ylabel('Electrode pair');
%     yticks(1:numel(cur_elec_contact_ind)-1);
%     set(gca, 'YTickLabel', pair_labels(yticks));
%     title(sprintf("PLV/ISPC adjacent contacts probe %s, %s", cur_letter, align_name));

%     nexttile;
%     delete(gca);
%     ax = gobjects(M,1);
%     top = 0.93; bottom = 0.09; left = 0.62; right = 0.94;
%     totalHeight = top - bottom;
%     h = totalHeight / M;  
%     for k = 1:M
%         % compute vertical position; want row 1 at top
%         pos = [left, top - k*h, right-left, h*0.95]; % spacing
%         ax(k) = axes('Position', pos);
%         plot(time, adj_pairs_lag_avg(k,:), 'k');
%         xlim([time(1) time(end)]);
%         ylim([min(adj_pairs_lag_avg(k,:)) max(adj_pairs_lag_avg(k,:))]); % use per-row limits
%         set(ax(k), 'YTick', []); 
%         if k < M
%             set(ax(k), 'XTick', [], 'XTickLabel', {});
%         end
%         % label each row on left
%         text(time(1) - 0.02*(time(end)-time(1)), mean(ylim), pair_labels(k), 'HorizontalAlignment','right', 'VerticalAlignment','middle');
%         box off;
%     end
%     %linkaxes(ax, 'x');
%     title(ax(1), 'Angle lag across time');
%     xlabel('Time (s)');
%     % check what y-lim defaults to (may need to set from -2pi to 2pi)

%     fname = sprintf("%s/probe%s_%s_avg_subplots.jpg", plot_out_dir, cur_letter, alignment);
%     print(fg1,'-djpeg',fname)
%     close;
% end

% ------------------------- 2. SINGLE-TRIAL ISPC -------------------------------

sesnum = 1;
for trial_num = 5:5:100
    ma_win = 50;

    for alignment = "fourth_unique_attribute" %alignments
        align_name = strrep(alignment, '_', ' ');
        adj_pairs_ISPC = []; % each row = time series of PLV for one pair of contacts
        adj_pairs_lag = []; % unwrap(angle1) - unwrap(angle2) at each timept; no binning/averaging across time
        pair_labels = strings(1, numel(cur_elec_contact_ind)-1);
        
        for cnum = 1:numel(cur_elec_contact_ind)-1 % for each adjacent pair of contacts
            [singleTrialISPC, angle_ts_lag] = ISPC_single_trial(subject_ID, reference, alignment, cur_letter, cnum, cnum+1, ma_win, trial_num, sesnum);
            adj_pairs_ISPC = [adj_pairs_ISPC; singleTrialISPC.'];
            adj_pairs_lag = [adj_pairs_lag; angle_ts_lag.'];
            pair_labels(cnum) = sprintf("%s-%s", cur_elec_contact_names(cnum), cur_elec_contact_names(cnum+1));
        end
        time = (1:width(adj_pairs_lag))/Fs;

        fg1 = figure('Units','normalized','Position',[0.1 0.1 0.8 0.5]); % subplots (left=ISPC as heatmap, right=lag as lineplot)
        tiledlayout(1,2,'TileSpacing','compact');
        
        nexttile;
        imagesc(time, 1:height(adj_pairs_ISPC), adj_pairs_ISPC); 
        colormap hot                   
        c = colorbar;
        c.Label.String = 'PLV'; 
        clim([0 1]);
        xlabel('Time (s)');
        ylabel('Electrode pair');
        yticks(1:numel(cur_elec_contact_ind)-1);
        set(gca, 'YTickLabel', pair_labels(yticks));
        title(sprintf("PLV/ISPC adjacent contacts probe %s, %s", cur_letter, align_name));

        nexttile;
        delete(gca);
        M = size(adj_pairs_lag,1); N = size(adj_pairs_lag,2);
        ax = gobjects(M,1);
        top = 0.93; bottom = 0.09; left = 0.62; right = 0.94;
        totalHeight = top - bottom;
        h = totalHeight / M;  
        for k = 1:M
            pos = [left, top - k*h, right-left, h*0.95];
            ax(k) = axes('Position', pos);
            plot(time, adj_pairs_lag(k,:), 'k');
            xlim([time(1) time(end)]);
            ylim([min(adj_pairs_lag(:)) max(adj_pairs_lag(:))]); % or use per-row limits [min(adj_pairs_lag(k,:)) max(adj_pairs_lag(k,:))]
            set(ax(k), 'YTick', []); 
            if k < M
                set(ax(k), 'XTick', [], 'XTickLabel', {});
            end
            text(time(1) - 0.02*(time(end)-time(1)), mean(ylim), pair_labels(k), 'HorizontalAlignment','right', 'VerticalAlignment','middle');
            box off;
        end
        title(ax(1), 'Angle lag across time');
        xlabel('Time (s)');

        final_out_dir = sprintf("%s/single_trial", plot_out_dir); if ~exist(final_out_dir, 'dir'), mkdir(final_out_dir); end
        fname = sprintf("%s/probe%s_%s_trial%d_subplots.jpg", final_out_dir, cur_letter, alignment, trial_num);
        print(fg1,'-djpeg',fname)
        close;
    end % end loop throguh alignments

end % end loop through trials


%% ISPC across trials
function [ISPC_avg, lag_avg] = ISPC_across_trials(subject_ID, reference, alignment, cur_letter, ct1, ct2)

% INPUTS: ct1, ct2 -- channel indices (into cur_elec_contact_ind)
% OUTPUT: Return time series of ISPC/PLV for specified alignment

addpath('/media/Data/Human_Intracranial_MAD/analysis/TravWaves/Code')
[~, data_base_dir, ~, ~, num_sessions, Fs, ~] = tw_setup(subject_ID, reference);
[cur_elec_contact_ind, ~] = get_single_probe_contacts(reference, subject_ID, cur_letter);

if strcmp(alignment, "trialstart") || strcmp(alignment, "anticipation")
    win = [Fs/2 0];
else
    win = [0 Fs/2];
end

inspect_alignments = ["first_unique_attribute", "second_unique_attribute", "third_unique_attribute", "fourth_unique_attribute", "inspection", "single_opt_first_inspection", "full_single_opt_info", "amount", "probability"];
if ismember(alignment, inspect_alignments)
    align_type = "inspection";
end

ISPC_ts = zeros(max(win), 1);
num_events = 0;
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

    % get xgp for ct1 and ct2
    contact1 = cur_elec_contact_ind(ct1);
    contact2 = cur_elec_contact_ind(ct2);
    if strcmp(reference,'Ground')
        load(sprintf('%s/%s/separate_channel_files/%s_MAD_SES%d_ch%03d.mat',data_base_dir,subject_ID,subject_ID,sesnum,contact1),"data"); 
        data1 = data;
        load(sprintf('%s/%s/separate_channel_files/%s_MAD_SES%d_ch%03d.mat',data_base_dir,subject_ID,subject_ID,sesnum,contact2),"data"); 
        data2 = data;
    elseif strcmp(reference,'neighbor_average')            
        load(sprintf('%s/%s/separate_channel_files/%s/%s_MAD_SES%d_ch%03d.mat',data_base_dir,subject_ID,reference,subject_ID,sesnum,contact1),"data");
        data1 = data;
        load(sprintf('%s/%s/separate_channel_files/%s/%s_MAD_SES%d_ch%03d.mat',data_base_dir,subject_ID,reference,subject_ID,sesnum,contact2),"data");
        data2 = data;
    end
        
    [b,a] = butter( 4, [5 40] ./ (Fs/2) ); filtered = filtfilt( b, a, data1 ); data1 = filtered;
    [b,a] = butter( 4, [5 40] ./ (Fs/2) ); filtered = filtfilt( b, a, data2 ); data2 = filtered;

    for event = 1:numel(align_times)
        data_padded1 = data1(round(align_times(event)-win(1)) - Fs:round(align_times(event)+win(2)-1) + Fs);
        data_padded2 = data2(round(align_times(event)-win(1)) - Fs:round(align_times(event)+win(2)-1) + Fs);
        
        xgp1 = generalized_phase_vector(data_padded1, Fs, 5); xgp1 = xgp1(1+Fs:1+Fs+win(2)-1);
        angle_ts1 = angle(xgp1); % time series of angles for ct1
        xgp2 = generalized_phase_vector(data_padded2, Fs, 5); xgp2 = xgp2(1+Fs:1+Fs+win(2)-1);
        angle_ts2 = angle(xgp2); % time series of angles for ct2

        ISPC_ts = ISPC_ts + exp((angle_ts1 - angle_ts2) * 1i);
    end
    num_events = num_events + numel(align_times);
end
ISPC_avg = abs(ISPC_ts / num_events);
lag_avg = angle(ISPC_ts / num_events);

end

%% ISPC for single trial
function [singleTrialISPC, angle_lag_ts] = ISPC_single_trial(subject_ID, reference, alignment, cur_letter, ct1, ct2, ma_win, trial_num, sesnum)

% INPUTS: ma_win (size of moving average window in samples, usually 50-150), trial_num
% OUTPUT: time series of ISPC/phase-locking (moving-averaged) for one trial
 
addpath('/media/Data/Human_Intracranial_MAD/analysis/TravWaves/Code')
[~, data_base_dir, ~, ~, ~, Fs, ~] = tw_setup(subject_ID, reference);
[cur_elec_contact_ind, ~] = get_single_probe_contacts(reference, subject_ID, cur_letter);

if strcmp(alignment, "trialstart") || strcmp(alignment, "anticipation")
    win = [Fs/2 0];
else
    win = [0 Fs/2];
end

inspect_alignments = ["first_unique_attribute", "second_unique_attribute", "third_unique_attribute", "fourth_unique_attribute", "inspection", "single_opt_first_inspection", "full_single_opt_info"];
if ismember(alignment, inspect_alignments)
    align_type = "inspection";
end

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

% get xgp for ct1 and ct2
contact1 = cur_elec_contact_ind(ct1);
contact2 = cur_elec_contact_ind(ct2);
if strcmp(reference,'Ground')
    load(sprintf('%s/%s/separate_channel_files/%s_MAD_SES%d_ch%03d.mat',data_base_dir,subject_ID,subject_ID,sesnum,contact1),"data"); 
    data1 = data;
    load(sprintf('%s/%s/separate_channel_files/%s_MAD_SES%d_ch%03d.mat',data_base_dir,subject_ID,subject_ID,sesnum,contact2),"data"); 
    data2 = data;
elseif strcmp(reference,'neighbor_average')            
    load(sprintf('%s/%s/separate_channel_files/%s/%s_MAD_SES%d_ch%03d.mat',data_base_dir,subject_ID,reference,subject_ID,sesnum,contact1),"data");
    data1 = data;
    load(sprintf('%s/%s/separate_channel_files/%s/%s_MAD_SES%d_ch%03d.mat',data_base_dir,subject_ID,reference,subject_ID,sesnum,contact2),"data");
    data2 = data;
end
    
[b,a] = butter( 4, [5 40] ./ (Fs/2) ); filtered = filtfilt( b, a, data1 ); data1 = filtered;
[b,a] = butter( 4, [5 40] ./ (Fs/2) ); filtered = filtfilt( b, a, data2 ); data2 = filtered;

data_padded1 = data1(round(align_times(trial_num)-win(1)) - Fs:round(align_times(trial_num)+win(2)-1) + Fs);
data_padded2 = data2(round(align_times(trial_num)-win(1)) - Fs:round(align_times(trial_num)+win(2)-1) + Fs);

xgp1 = generalized_phase_vector(data_padded1, Fs, 5); xgp1 = xgp1(1+Fs:1+Fs+win(2)-1);
angle_ts1 = angle(xgp1); % time series of angles for ct1
xgp2 = generalized_phase_vector(data_padded2, Fs, 5); xgp2 = xgp2(1+Fs:1+Fs+win(2)-1);
angle_ts2 = angle(xgp2); % time series of angles for ct2

complex_ts = movmean(exp((angle_ts1 - angle_ts2) * 1i), ma_win);
angle_lag_ts = angle(exp(1i * (angle_ts1 - angle_ts2)));
singleTrialISPC = abs(complex_ts);


% switch freq_name
%    case 'Delta'
%        binSize = 100; % >10 samples per delta cycle
%    case 'Theta'
%        binSize = 30; % 10 samples per theta cycle
%    case 'Alpha'
%        binSize = 20; % 10 samples per alpha cycle 
%    case 'Beta'
%        binSize = 10; % 9-10 samples per beta cycle
%    case 'Gamma_joint'
%        binSize = 4;
%    case 'Gamma_high'
%        binSize = 4; % 4-5 samples per high gamma cycle
%    case 'Gamma_low'
%        binSize = 7; % 6 samples per low gamma cycle
% end

end % end function

%% Phase lag/latency
% function phase_lag_ts = phase_lag(signal1, signal2, Fs)
%     % INPUT: real-valued signal -- clipped and filtered -- for 2 channels
%     % OUTPUT: mean phase lag in radians (and milliseconds)
%     % - Positive lag = signal2 lags signal1 (wave travels 1→2)

%     analytic1 = generalized_phase_vector(signal1, Fs, 5);
%     analytic2 = generalized_phase_vector(signal2, Fs, 5);
    
%     phase_lag_ts = angle(analytic1) - angle(analytic2);
    
%     % Circular mean
%     %mean_phase_lag = angle(mean(exp(1i * phase_diff), 1));
    
% end