% run_phase_corr.m -- 3/15/26
% Use function phase_dist_corr to compute circular-linear corr between contact position on probe and phase angle
% Plot time series of phase-dist correlations

% Check if there are large chunks of gray/white matter; if so, run phase_dist_corr on those chunks
% 1. phase-dist corr on chunks of wm/gm + everything else
% 2. phase-dist corr on entire probes
% 3. phase-dist corr on medial half and lateral half of probe

% "parameters"
subject_ID = 'EMU038';
reference = 'neighbor_average';
cur_letter = "N";
average_across_trials = true;

addpath('/media/Data/Human_Intracranial_MAD/analysis/TravWaves/Code')
addpath('/media/Data/Human_Intracranial_MAD/analysis/TravWaves/Code/ERP')
addpath('/media/Data/Human_Intracranial_MAD/analysis/PAC_code/CF_Coupling/generalized phase')
[datapath, data_base_dir, tw_out_dir, out_dir, num_sessions, Fs, elec_letters] = tw_setup(subject_ID, reference);
[cur_elec_contact_ind, cur_elec_contact_names] = get_single_probe_contacts(reference, subject_ID, cur_letter);
[wm_gm_chunks, chunk_names, chunk_areas] = gm_wm_chunking(subject_ID, reference, cur_letter);

win = [0, Fs/2];
subwin_st = win(1) + 1; % compute corr over entire event window
subwin_end = win(2);

if average_across_trials
    save_dir_base = sprintf("%s/Plots/%s/%s/circ_lin_corr/Probe %s/trial_avg", tw_out_dir, reference, subject_ID, cur_letter);
else
    save_dir_base = sprintf("%s/Plots/%s/%s/circ_lin_corr/Probe %s", tw_out_dir, reference, subject_ID, cur_letter);
end
if ~exist(save_dir_base, 'dir'), mkdir(save_dir_base); end

%alignment = "third_unique_attribute"; % "fourth_unique_attribute"

if average_across_trials
    fg1 = figure;
    set( fg1, 'position', [ 88  100  1450  760 ] )
end
for alignment = ["first_unique_attribute", "second_unique_attribute","third_unique_attribute","fourth_unique_attribute"] %["full_single_opt_info","single_opt_first_inspection"] %

    align_names = strrep(alignment, '_', ' '); align_type = alignment;
    [~, avg_events_allct1, ~, ~, ~, ~, ~, ~] = plot_save_trial_avg(subject_ID, reference, align_type, alignment, align_names, cur_letter, false, false, false);
    
    % ------------ Specific probe subsets ----------------
    % based on large consecutive chunks of gray/white matter

    % if any(diff(wm_gm_chunks) > 3)
    %     save_dir = save_dir_base;
    %     large_chunks_st = wm_gm_chunks(diff(wm_gm_chunks) > 3);
    %     large_chunks_names = chunk_names(diff(wm_gm_chunks) > 3);
    %     large_chunks_areas = chunk_areas(diff(wm_gm_chunks) > 3);
    %     for chunk = 1:numel(large_chunks_st)
    %         idx = find(wm_gm_chunks == large_chunks_st(chunk));
    %         ch_nums = large_chunks_st(chunk):wm_gm_chunks(idx+1)-1;

    %         % ch_nums = 7:12; % TEMP manual
            
    %         allses_angle_cts = compute_angle_ts(subject_ID, reference, cur_letter, ch_nums, alignment, subwin_st, subwin_end);
    %         [cur_chunk_corr, ~] = phase_dist_corr(ch_nums, subwin_st, subwin_end, allses_angle_cts);
    %         % cur_chunk_pv -- to use this var or not?

    %         % Plot time series of correlations (cur_chunk_corr)
    %             % plot non-significant corr-coefs in blue, significant corrs in red 
    %             % ... but need multiple comparisons correction!

    %         area = large_chunks_areas(chunk); % WM/GM/other
    %         name_st = large_chunks_names(chunk);
    %         idx_name = find(strcmp(chunk_names, name_st));
    %         name_end = chunk_names(idx_name+1);

    %         % name_st = "N6"; name_end = "N1"; % TEMP manual

    %         if average_across_trials
    %             avg_corr_ts = mean(cur_chunk_corr, 2);
                
    %             k = 50; 
    %             p1 = plot((subwin_st:subwin_end)/Fs, movmean(avg_corr_ts, k),'LineWidth',2, 'DisplayName', sprintf("%s smoothed", align_names));
    %             hold on
    %             plot((subwin_st:subwin_end)/Fs, avg_corr_ts, 'Color', [p1.Color 0.3],'LineWidth',2, 'DisplayName', sprintf("%s raw", align_names));
    %             xlabel('time from alignment (s)'); 
    %             ylabel('corr'); %ylim(y1);
    %         else
    %             for trial = 1:5:100 %width(cur_chunk_corr) % TEMP 3/13/26
    %                 fg1 = figure;
    %                 plot((subwin_st:subwin_end)/Fs, cur_chunk_corr(:,trial))
    %                 title(sprintf("Circular-Linear Correlation -- %s %s-%s %s", align_names, name_st, name_end, area))
    %                 xlabel('time from alignment (s)'); 
    %                 ylabel('corr');
    %                 ylim([0 1.0]);
    %                 fname = sprintf("%s/%s_%s-%s_%s_trial%d.jpg",save_dir_base, alignment, name_st, name_end, area, trial);
    %                 print(fg1,'-djpeg',fname)
    %                 close;
    %             end
    %         end % end if average_across_trials

    %     end
        
    % end


    % ----------- Generic probe subsets ------------
    probe_subsets = ["whole","medial","lateral"];
    probe_subset = "whole";

    if strcmp(probe_subset, "whole")
        ch_nums = 1:numel(cur_elec_contact_ind);
        name_st = cur_elec_contact_names(1);
        name_end = cur_elec_contact_names(end);
        save_dir = sprintf("%s/whole_probe", save_dir_base);
    elseif strcmp(probe_subset, "medial")
        ch_nums = 1:ceil((numel(cur_elec_contact_ind))/2);
        name_st = cur_elec_contact_names(end);
        name_end = cur_elec_contact_names(end - ceil((numel(cur_elec_contact_ind))/2) + 1);
        save_dir = sprintf("%s/medial_half_probe_cts", save_dir_base);
    elseif strcmp(probe_subset, "lateral")
        ch_nums = ceil((numel(cur_elec_contact_ind))/2)+1:numel(cur_elec_contact_ind);
        name_st = cur_elec_contact_names(end - ceil((numel(cur_elec_contact_ind))/2));
        name_end = cur_elec_contact_names(1);
        save_dir = sprintf("%s/lateral_half_probe_cts", save_dir_base);
    end
    if ~exist(save_dir, 'dir'), mkdir(save_dir); end

    allses_angle_cts = compute_angle_ts(subject_ID, reference, cur_letter, ch_nums, alignment, subwin_st, subwin_end);
    [cur_chunk_corr, cur_chunk_pv] = phase_dist_corr(ch_nums, subwin_st, subwin_end, allses_angle_cts);

    area = probe_subset;
    if average_across_trials
        avg_corr_ts = mean(cur_chunk_corr, 2);

        k = 50; 
        p1 = plot((subwin_st:subwin_end)/Fs, movmean(avg_corr_ts, k),'LineWidth',2, 'DisplayName', sprintf("%s smoothed", align_names));
        hold on
        plot((subwin_st:subwin_end)/Fs, avg_corr_ts, 'Color', [p1.Color 0.3],'LineWidth',2, 'DisplayName', sprintf("%s raw", align_names));
        xlabel('time from alignment (s)'); 
        ylabel('corr');

    else
        for trial = 1:5:100 %width(cur_chunk_corr)
            fg1 = figure;
            set( fg1, 'position', [ 88  100  1450  760 ] )
            tiledlayout(1,2,'TileSpacing','compact');
            title(sprintf("Circular-Linear Correlation -- %s %s-%s %s", align_names, name_st, name_end, area))

            % Raw -- full time resolution
            ax1 = nexttile;
            plot((subwin_st:subwin_end)/Fs, cur_chunk_corr(:,trial)) % '.', 'MarkerSize', 8
            title(ax1, 'Raw (full time resolution)')
            xlabel('time from alignment (s)'); 
            ylabel('corr');
            ylim([0 1.0]);

            % Moving average -- lower/smoothed time resolution
            ax2 = nexttile;
            k = 50;
            plot((subwin_st:subwin_end)/Fs, movmean(cur_chunk_corr(:,trial), k))
            title(ax2, sprintf("Moving average (k = %d smoothed)", k))
            xlabel('time from alignment (s)'); 
            ylabel('corr');
            ylim([0 1.0]);

            fname = sprintf("%s/%s_%s-%s_%s_trial%d.jpg",save_dir, alignment, name_st, name_end, area, trial);
            print(fg1,'-djpeg',fname)
            close;
        end

    end
end % end loop thru alignments
if average_across_trials
    legend;
    title(sprintf("Circular-Linear Correlation Time Series"))
    fname = sprintf("%s/all4inspect_%s-%s_%s.jpg",save_dir, name_st, name_end, area); %full_vs_partial_info
    print(fg1,'-djpeg',fname)
    close;
end