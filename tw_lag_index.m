% tw_lag_index.m: plot traveling wave index (proportion of time points with spatiotemporally ordered lags) over time
addpath('/media/Data/Human_Intracranial_MAD/analysis/TravWaves/Code')
subject_ID = 'EMU038'; reference = 'neighbor_average';
ma_win = 0.01; cur_letter = 'N'; chnums = 8:10; %10:-1:8
subwin_st = 150; subwin_end = 300;

[~, ~, ~, ~, ~, Fs, ~] = tw_setup(subject_ID, reference);
[cur_elec_contact_ind, cur_elec_contact_names] = get_single_probe_contacts(reference, subject_ID, cur_letter);
ct1 = cur_elec_contact_names(chnums(end)); ct2 = cur_elec_contact_names(chnums(1));

alignment = "fourth_unique_attribute"; %fourth_unique_attribute

lag_index_ts = phase_lag_index(subject_ID, reference, cur_letter, chnums, ma_win, alignment, subwin_st, subwin_end);
avg_lag_index = mean(lag_index_ts, 1);

% save_dir = sprintf('/media/Data/Human_Intracranial_MAD/analysis/TravWaves/Results/%s/%s/Probe %s', reference, subject_ID, cur_letter);
% fname = sprintf('%s/%s_%s-%s_phase_in_order_bool.mat', save_dir, alignment, ct1, ct2);
% save(fname, "lag_index_ts")

% % ----- HEATMAP HISTOGRAMS -----
% % each column is a histogram showing distribution of proportions across trials for a given timepoint in the event window
% numBins = 20;
% binEdges = linspace(0, 1, numBins + 1);
% binCenters = binEdges(1:end-1) + diff(binEdges)/2;

% heatmapData = zeros(numBins, size(lag_index_ts, 2));
% for i = 1:size(lag_index_ts, 2)
%     heatmapData(:, i) = histcounts(lag_index_ts(:, i), binEdges, "Normalization","probability");
% end

% figure;
% imagesc((1:size(lag_index_ts, 2))/Fs, binCenters, heatmapData);
% set(gca, 'YDir', 'normal'); % flip y-axis so 0 is at the bottom
% colormap(parula);           % 'hot', 'magma'
% colorbar;

% xlabel('Timepoint (s)');
% ylabel('Phase Lag Proportion');
% title(sprintf("Phase Lag Index Distribution Across Trials %s-%s", ct1, ct2));


% % ----- LINE PLOT -----
% figure; 
% % plot((1:numel(avg_lag_index))/Fs, avg_lag_index, 'LineWidth', 2)
% plot((1:numel(avg_lag_index))/Fs, lag_index_ts(1,:), 'LineWidth', 1.5)
% hold on;
% for trial = 5:10
%     plot((1:numel(avg_lag_index))/Fs, lag_index_ts(trial,:), 'LineWidth', 1.5)
% end
% xlabel("time from alignment (s)")
% ylabel("proportion ordered phase lag")
% %ylim([0, 1])
% title(sprintf("Phase Lag Index %s-%s", ct1, ct2))
% %close;

% --- FUNCTION DEFINITION ----
function lag_index_ts = phase_lag_index(subject_ID, reference, cur_letter, chnums, ma_win, alignment, subwin_st, subwin_end)

% INPUTS
% - chnums: channel numbers (indices into cur_elec_contact_ind) as an int array, in the order in which wave should be detected
%           e.g., 5:7 --> check if ch5 < ch6 < ch7 (5 lags behind 6, 6 lags behind 7)
% - ma_win: width of time window (in seconds) for computing proportion/count of ordered lags
% - subwin_st / end: start and end times (in samples/indices) for the window on which to compute the index

% OUTPUT: N_trials x N_timepts matrix, each row is time series of proportions of timepts in each moving window
%         at which contacts lag behind one another in phase in the correct sequence, NOT trial-averaged

addpath('/media/Data/Human_Intracranial_MAD/analysis/TravWaves/Code')
addpath('/media/Data/Human_Intracranial_MAD/analysis/TravWaves/Code/ERP')
addpath('/media/Data/Human_Intracranial_MAD/analysis/TravWaves/Code/generalized-phase/analysis')
[~, data_base_dir, ~, ~, num_sessions, Fs, ~] = tw_setup(subject_ID, reference);
[cur_elec_contact_ind, ~] = get_single_probe_contacts(reference, subject_ID, cur_letter);
ma_win = round(ma_win * Fs);

if strcmp(alignment, "trialstart") || strcmp(alignment, "anticipation")
    win = [Fs/2 0];
else
    win = [0 Fs/2];
end

inspect_alignments = ["first_unique_attribute", "second_unique_attribute", "third_unique_attribute", "fourth_unique_attribute", "inspection", "single_opt_first_inspection", "full_single_opt_info"];
if ismember(alignment, inspect_alignments)
    align_type = "inspection";
end

lag_index_ts = []; % each row is a time series of proportions

% create binary/boolean array of 1s/0s, then apply moving avg of size win_size on that array -> will return a proportion
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

    all_ch_filt = []; % each col = 1 channel's time series
    for cnum = chnums
        contact = cur_elec_contact_ind(cnum);
        if strcmp(reference,'Ground')
            load(sprintf('%s/%s/separate_channel_files/%s_MAD_SES%d_ch%03d.mat',data_base_dir,subject_ID,subject_ID,sesnum,contact),"data"); 
        elseif strcmp(reference,'neighbor_average')            
            load(sprintf('%s/%s/separate_channel_files/%s/%s_MAD_SES%d_ch%03d.mat',data_base_dir,subject_ID,reference,subject_ID,sesnum,contact),"data");
        end
        [b,a] = butter( 4, [5 40] ./ (Fs/2) ); filtered = filtfilt( b, a, data ); data = filtered;
        all_ch_filt = [all_ch_filt data];
    end

    for event = 1:numel(align_times)
        
        all_ch_event_padded = all_ch_filt(round(align_times(event)-win(1)) - Fs:round(align_times(event)+win(2)-1) + Fs, :);
        %unpadded = all_ch_event_padded(1+Fs:1+Fs+win(2)-1, :);
        %all_ch_event_zscored = (unpadded - mean(unpadded, 1)) ./ std(unpadded, 0, 1);

        all_ch_angle_event = [];
        for cnum_idx = 1:numel(chnums)
            % apply xgp to each column in a padded event windopw
            xgp_padded = generalized_phase_vector(all_ch_event_padded(:,cnum_idx), Fs, 5);
            xgp_event = xgp_padded(1+Fs:1+Fs+win(2)-1);
            angle_event = angle(xgp_event);
            all_ch_angle_event = [all_ch_angle_event angle_event(subwin_st:subwin_end)];
        end

        % % REFERENCE time point:
        % [~, linearIdx] = max(all_ch_event_padded, [], 'all');
        % [t0, ~] = ind2sub(size(A), linearIdx); % t0 is reference time

        % ------------------
        % % METHOD 1 -- don't unwrap angles
        % wrapped_angle_diffs = zeros(height(all_ch_angle_event), numel(chnums)-1);
        % for cnum_idx = 1:numel(chnums)-1
        %     angle_diff = all_ch_angle_event(:,cnum_idx+1) - all_ch_angle_event(:,cnum_idx);
        %     %wrapped_angle_diffs(:,cnum) = atan(sin(angle_diff), cos(angle_diff));
        %     for tp = 1:height(wrapped_angle_diffs)
        %         wrapped_angle_diffs(tp, cnum) = atan(sin(angle_diff(tp)) / cos(angle_diff(tp)));
        %     end
        % end
        % bool_ts_event = (wrapped_angle_diffs.' > 0); % transpose so that time is along horizontal

        % METHOD 2 -- unwrap angles
        unwrapped_angle_diffs = zeros(height(all_ch_angle_event), numel(chnums)-1);
        for cnum = 1:numel(chnums)-1
            % angle_diff = all_ch_angle_event(:,cnum+1) - all_ch_angle_event(:,cnum);
            % unwrapped_angle_diffs(:,cnum) = unwrap(angle(exp(1i * angle_diff)));

            unwrapped_angle_diffs(:,cnum) = unwrap(all_ch_angle_event(:,cnum+1)) - unwrap(all_ch_angle_event(:,cnum));

            % for pc = 1:100:height(unwrapped_angle_diffs)
            %     if pc+100 > height(unwrapped_angle_diffs)
            %         unwrapped_angle_diffs(pc:end, cnum) = unwrap(all_ch_angle_event(pc:end,cnum+1)) - unwrap(all_ch_angle_event(pc:end,cnum));
            %     else
            %         unwrapped_angle_diffs(pc:pc+100, cnum) = unwrap(all_ch_angle_event(pc:pc+100,cnum+1)) - unwrap(all_ch_angle_event(pc:pc+100,cnum));
            %     end
            % end
        end
        bool_ts_event = (unwrapped_angle_diffs.' > 0);
        % -----------------
        
        bool_ts_event = sum(bool_ts_event, 1) == numel(chnums)-1; % all chnums are in order
            % nnz(bool_ts_event) is number of time points at which phase of chnums are in order
        %prop_lag_order = movmean(bool_ts_event, ma_win, 2);
        %lag_index_ts = [lag_index_ts; prop_lag_order];
        lag_index_ts = [lag_index_ts; bool_ts_event];

    end % end loop through events

end


end