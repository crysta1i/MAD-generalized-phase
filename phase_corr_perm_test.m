% phase_corr_perm_test.m
% Test significance of circular-linear corr coefs (section 1) and differences of corr coefs (section 2)
% 1000 random permutations of the contacts on the probe -- apply each to all timepoints in each trial

%% 0. Setup -- always run this section

subject_ID = 'EMU038';
reference = 'neighbor_average';
cur_letter = "M";
average_across_trials = true;
iterations = 1000;
alignments = ["single_opt_first_inspection","full_single_opt_info"]; 
align_names = ["Single Option First Inspection","Full Single Option Info"];

addpath('/media/Data/Human_Intracranial_MAD/analysis/TravWaves/Code')
addpath('/media/Data/Human_Intracranial_MAD/analysis/TravWaves/Code/ERP')
addpath('/media/Data/Human_Intracranial_MAD/analysis/PAC_code/CF_Coupling/generalized phase')
[datapath, data_base_dir, tw_out_dir, out_dir, num_sessions, Fs, elec_letters] = tw_setup(subject_ID, reference);
[cur_elec_contact_ind, cur_elec_contact_names] = get_single_probe_contacts(reference, subject_ID, cur_letter);
[wm_gm_chunks, chunk_names, chunk_areas] = gm_wm_chunking(subject_ID, reference, cur_letter);

win = [0, Fs/2];
subwin_st = win(1) + 1; % compute corr over entire event window
subwin_end = win(2);
ch_nums = 1:numel(cur_elec_contact_ind); name_st = cur_elec_contact_names(end); name_end = cur_elec_contact_names(1);

if average_across_trials
    save_dir_base = sprintf("%s/Plots/%s/%s/circ_lin_corr/Probe %s/trial_avg", tw_out_dir, reference, subject_ID, cur_letter);
else
    %save_dir_base = sprintf("%s/Plots/%s/%s/circ_lin_corr/Probe %s", tw_out_dir, reference, subject_ID, cur_letter);
end
if ~exist(save_dir_base, 'dir'), mkdir(save_dir_base); end
save_dir = sprintf("%s/whole_probe", save_dir_base);

%% 1. Testing significance of corr at each timepoint for each alignment separately
fg1 = figure;
set( fg1, 'position', [ 88  100  1450  760 ] )

for align_num = 1:numel(alignments)
    alignment = alignments(align_num);
    rng(10); % does this apply same permutations across the alignments?
    
    %if average_across_trials
        
        corr_perm_distr = zeros(iterations, abs(win(2) - win(1))); % one permutation distr per timept - each column=distr for 1 timept
        allses_angle_cts = compute_angle_ts(subject_ID, reference, cur_letter, ch_nums, alignment, subwin_st, subwin_end);
        [obs_corr, ~] =  phase_dist_corr(ch_nums, subwin_st, subwin_end, allses_angle_cts);
        obs_corr = mean(obs_corr, 2);
        for iter = 1:iterations
            allses_angle_cts_perm = allses_angle_cts(:, randperm(numel(ch_nums)), :); % applies same permutation across all trials
            [corr_mat, ~] = phase_dist_corr(ch_nums, subwin_st, subwin_end, allses_angle_cts_perm);

            avg_corr_ts_perm = mean(corr_mat, 2);
            corr_perm_distr(iter, :) = avg_corr_ts_perm;
            
        end

        sig = zeros(1, width(corr_perm_distr)); % boolean time series indicating which timepoints have significant corr
        for tp = 1:width(corr_perm_distr)
            cur_tp_distr = sort(corr_perm_distr(:,tp),"ascend");
            threshold = round(0.95 * iterations); % one-sided test
            if obs_corr(tp) > cur_tp_distr(threshold)
                sig(tp) = 1;
            end
        end
    %else
        % for trial = 1:5:100
        %     % ch_nums_perm = randperm(ch_nums);
        %     % applies different permutation for each trial
        % end
    %end
    k=50;
    p1 = plot((subwin_st:subwin_end)/Fs, movmean(obs_corr, k),'LineWidth',2, 'DisplayName', sprintf("%s smoothed", align_names(align_num)));
    hold on
    plot((subwin_st:subwin_end)/Fs, obs_corr, 'Color', [p1.Color 0.3],'LineWidth',2, 'DisplayName', sprintf("%s raw", align_names(align_num)));
    xlabel('time from alignment (s)'); 
    ylabel('corr');
    y_lims = ylim;

    % ------ shade regions of plot that are stat significant -------
        % look for consecutive runs in sig: 
    start_sig_idx = find(diff(sig) > 0); end_sig_idx = find(diff(sig) < 0);
    if isempty(start_sig_idx), start_sig_idx = 1; end
    if isempty(end_sig_idx), end_sig_idx = numel(sig); end
    if start_sig_idx(1) > end_sig_idx(1), start_sig_idx = [1 start_sig_idx]; end
    if start_sig_idx(end) > end_sig_idx(end), end_sig_idx = [end_sig_idx numel(sig)]; end
    for idx = 1:numel(start_sig_idx)
        x_start = (start_sig_idx(idx)-1) / Fs; % Start time
        x_end = end_sig_idx(idx) / Fs;   % End time

        % Define the 4 corners of the rectangle
        x_coords = [x_start, x_end, x_end, x_start];
        y_coords = [y_lims(1), y_lims(1), y_lims(2), y_lims(2)];

        % Draw the shaded region
        p = patch(x_coords, y_coords, [p1.Color]); 
        p.FaceAlpha = 0.1;  % Set transparency (0 to 1)
        p.EdgeColor = 'none';

        uistack(p, 'bottom'); % Ensure the data plot stays on top
    end
    % -------------------------------------------------------------

end
legend;
title(sprintf("Circular-Linear Correlation Time Series"))
fname = sprintf("%s/full_vs_partial_info_%s-%s_show_sig.jpg",save_dir, name_st, name_end); %full_vs_partial_info % all4inspect
print(fg1,'-djpeg',fname)
close;

%% 2. Testing significnace of DIFFERENCE in corr between alignments (single_opt_first_inspection and full_single_opt_info)
% create permutation distr of corr_1st_inspection - corr_2nd_inspection

corr_diff_perm_distr = zeros(iterations, abs(win(2) - win(1))); 
sig_diff = zeros(1, width(corr_diff_perm_distr));

rng(10);

allses_angle_cts1 = compute_angle_ts(subject_ID, reference, cur_letter, ch_nums, "single_opt_first_inspection", subwin_st, subwin_end);
allses_angle_cts2 = compute_angle_ts(subject_ID, reference, cur_letter, ch_nums, "full_single_opt_info", subwin_st, subwin_end);
[obs_corr1, ~] = phase_dist_corr(ch_nums, subwin_st, subwin_end, allses_angle_cts1);
[obs_corr2, ~] = phase_dist_corr(ch_nums, subwin_st, subwin_end, allses_angle_cts2);
obs_diff = obs_corr1 - obs_corr2;
obs_diff = mean(obs_diff, 2);
for iter = 1:iterations
    ch_nums_perm = randperm(numel(ch_nums));
    allses_angle_cts_perm1 = allses_angle_cts1(:, ch_nums_perm, :);
    allses_angle_cts_perm2 = allses_angle_cts2(:, ch_nums_perm, :);

    [cur_perm_corr, ~] = phase_dist_corr(ch_nums, subwin_st, subwin_end, allses_angle_cts_perm1);
    avg_corr_ts_perm1 = mean(cur_perm_corr, 2);
    
    [cur_perm_corr, ~] = phase_dist_corr(ch_nums, subwin_st, subwin_end, allses_angle_cts_perm2);
    avg_corr_ts_perm2 = mean(cur_perm_corr, 2);

    corr_diff_perm_distr(iter, :) = avg_corr_ts_perm1 - avg_corr_ts_perm2;
end

for tp = 1:width(corr_diff_perm_distr)
    cur_diff_distr = sort(corr_diff_perm_distr(:,tp), "ascend");
    threshold_upper = round(0.975 * iterations); % two-sided test
    threshold_lower = round(0.025 * iterations);
    if obs_diff(tp) > cur_diff_distr(threshold_upper) || obs_diff(tp) < cur_diff_distr(threshold_lower)
        sig_diff(tp) = 1;
    end
end

fg1 = figure;
set( fg1, 'position', [ 88  100  1450  760 ] )

k=50;
p1 = plot((subwin_st:subwin_end)/Fs, movmean(mean(obs_corr1,2), k),'LineWidth',2, 'DisplayName', "Single opt First Inspection smoothed");
hold on
plot((subwin_st:subwin_end)/Fs, mean(obs_corr1, 2), 'Color', [p1.Color 0.3],'LineWidth',2, 'DisplayName', "Single opt First Inspection raw");
p2 = plot((subwin_st:subwin_end)/Fs, movmean(mean(obs_corr2, 2), k),'LineWidth',2, 'DisplayName', "Full Single Opt Info smoothed");
plot((subwin_st:subwin_end)/Fs, mean(obs_corr2, 2), 'Color', [p2.Color 0.3],'LineWidth',2, 'DisplayName', "Full Single Opt Info raw");
xlabel('time from alignment (s)'); 
ylabel('corr');
y_lims = ylim;

start_sig_idx = find(diff(sig_diff) > 0); end_sig_idx = find(diff(sig_diff) < 0);
if isempty(start_sig_idx), start_sig_idx = 1; end
if isempty(end_sig_idx), end_sig_idx = numel(sig_diff); end
if start_sig_idx(1) > end_sig_idx(1), start_sig_idx = [1 start_sig_idx]; end
if start_sig_idx(end) > end_sig_idx(end), end_sig_idx = [end_sig_idx numel(sig_diff)]; end
for idx = 1:numel(start_sig_idx)
    x_start = (start_sig_idx(idx)-1) / Fs; % Start time
    x_end = end_sig_idx(idx) / Fs;   % End time

    % Define the 4 corners of the rectangle
    x_coords = [x_start, x_end, x_end, x_start];
    y_coords = [y_lims(1), y_lims(1), y_lims(2), y_lims(2)];

    % Draw the shaded region
    p = patch(x_coords, y_coords, [p2.Color]); % 'DisplayName', "Sig Difference in Corr"
    p.FaceAlpha = 0.1;  % Set transparency (0 to 1)
    p.EdgeColor = 'none';

    uistack(p, 'bottom'); % Ensure the data plot stays on top
    if idx == 1
        p.DisplayName = "Sig Difference in Corr";
    else
        p.HandleVisibility = 'off';
    end
end

legend;
title(sprintf("Circular-Linear Correlation Time Series"))
fname = sprintf("%s/full_vs_partial_info_%s-%s_show_sig.jpg",save_dir, name_st, name_end); %full_vs_partial_info %all4inspect
print(fg1,'-djpeg',fname)
close;


%% 3. Testing for significant (paired) difference between 1st and 2nd inspection of option correlation
% Saves results (p-values and clusters) of cluster-based permutation test (using permutest.m function)

addpath('/media/Data/Human_Intracranial_MAD/_toolbox/spectral-analysis/cluster_statistic_code/')

save_dir = sprintf('/media/Data/Human_Intracranial_MAD/analysis/TravWaves/Results/%s/%s/Probe %s', reference, subject_ID, cur_letter);
if ~exist(save_dir, 'dir'), mkdir(save_dir); end

allses_angle_cts = compute_angle_ts(subject_ID, reference, cur_letter, ch_nums, 'single_opt_first_inspection', subwin_st, subwin_end);
[first_inspection_correlations, ~] = phase_dist_corr(ch_nums, subwin_st, subwin_end, allses_angle_cts);
allses_angle_cts = compute_angle_ts(subject_ID, reference, cur_letter, ch_nums, 'full_single_opt_info', subwin_st, subwin_end);
[second_inspection_correlations, ~] = phase_dist_corr(ch_nums, subwin_st, subwin_end, allses_angle_cts);
if width(first_inspection_correlations) > width(second_inspection_correlations)
    first_inspection_correlations = first_inspection_correlations(:,1:width(second_inspection_correlations)); 
elseif width(first_inspection_correlations) < width(second_inspection_correlations)
    second_inspection_correlations = second_inspection_correlations(:,1:width(first_inspection_correlations)); 
end
[clusters, p_values, t_sums, permutation_distribution] = permutest(first_inspection_correlations, second_inspection_correlations, 1, 0.05, iterations, 1);

fname = sprintf('%s/first_second_inspect_cluster_perm_results.mat', save_dir);
save(fname, "clusters", "p_values", "t_sums", "permutation_distribution")

%% 4. Plotting based on significant time points from cluster perm test

base_dir = sprintf('/media/Data/Human_Intracranial_MAD/analysis/TravWaves/Results/%s/%s/Probe %s', reference, subject_ID, cur_letter);
result_path = sprintf('%s/first_second_inspect_cluster_perm_results.mat', base_dir);
load(result_path, "clusters", "p_values")

allses_angle_cts1 = compute_angle_ts(subject_ID, reference, cur_letter, ch_nums, "single_opt_first_inspection", subwin_st, subwin_end);
allses_angle_cts2 = compute_angle_ts(subject_ID, reference, cur_letter, ch_nums, "full_single_opt_info", subwin_st, subwin_end);
[obs_corr1, ~] = phase_dist_corr(ch_nums, subwin_st, subwin_end, allses_angle_cts1);
[obs_corr2, ~] = phase_dist_corr(ch_nums, subwin_st, subwin_end, allses_angle_cts2);

start_sig_idx = []; end_sig_idx = [];
for i = 1:numel(p_values)
    if p_values(i) < 0.05
        sig_time_idx = clusters{i};
        start_sig_idx = [start_sig_idx; sig_time_idx(1)];
        end_sig_idx = [end_sig_idx; sig_time_idx(end)];
    end
end

fg1 = figure;
set( fg1, 'position', [ 88  100  1450  760 ] )

k=50;
p1 = plot((subwin_st:subwin_end)/Fs, movmean(mean(obs_corr1,2), k),'LineWidth',2, 'DisplayName', "Single opt First Inspection smoothed");
hold on
plot((subwin_st:subwin_end)/Fs, mean(obs_corr1, 2), 'Color', [p1.Color 0.3],'LineWidth',2, 'DisplayName', "Single opt First Inspection raw");
p2 = plot((subwin_st:subwin_end)/Fs, movmean(mean(obs_corr2, 2), k),'LineWidth',2, 'DisplayName', "Full Single Opt Info smoothed");
plot((subwin_st:subwin_end)/Fs, mean(obs_corr2, 2), 'Color', [p2.Color 0.3],'LineWidth',2, 'DisplayName', "Full Single Opt Info raw");
xlabel('time from alignment (s)'); 
ylabel('corr');
y_lims = ylim;

if ~isempty(start_sig_idx)
    for idx = 1:numel(start_sig_idx)
        x_start = start_sig_idx(idx) / Fs; % Start time
        x_end = end_sig_idx(idx) / Fs;   % End time

        % Define the 4 corners of the rectangle
        x_coords = [x_start, x_end, x_end, x_start];
        y_coords = [y_lims(1), y_lims(1), y_lims(2), y_lims(2)];

        % Draw the shaded region
        p = patch(x_coords, y_coords, [p2.Color]); % 'DisplayName', "Sig Difference in Corr"
        p.FaceAlpha = 0.1;  % Set transparency (0 to 1)
        p.EdgeColor = 'none';

        uistack(p, 'bottom'); % Ensure the data plot stays on top
        if idx == 1
            p.DisplayName = "Sig Difference in Corr";
        else
            p.HandleVisibility = 'off';
        end
    end
end

legend;
title(sprintf("Circular-Linear Correlation Time Series"))
fname = sprintf("%s/full_vs_partial_info_%s-%s_cluster_permtest_sig.jpg",save_dir, name_st, name_end); %full_vs_partial_info %all4inspect
print(fg1,'-djpeg',fname)
close;