% Inter-trial phase clustering analysis - 3/17/26
% Check how consistent phase angles are at various time points in a window
% by computing and plotting a time series of ITPC values (avg across trials)
% Line 81 -- ITPC function def'n 
% --------------------------------------------

addpath('/media/Data/Human_Intracranial_MAD/analysis/TravWaves/Code')
subject_ID = 'EMU001'; reference = 'neighbor_average';
alignments = ["amount","probability"]; %["first_unique_attribute", "second_unique_attribute","third_unique_attribute","fourth_unique_attribute"];

[~, ~, tw_out_dir, ~, num_sessions, Fs, elec_letters] = tw_setup(subject_ID, reference);

for cur_letter = elec_letters.'
    [cur_elec_contact_ind, cur_elec_contact_names] = get_single_probe_contacts(reference, subject_ID, cur_letter);
    plot_out_dir = sprintf("%s/Plots/%s/%s/ITPC/Probe %s", tw_out_dir, reference, subject_ID, cur_letter);
    if ~exist(plot_out_dir, 'dir'), mkdir(plot_out_dir); end

    if numel(cur_elec_contact_ind) == 18
        cnum_runs = {1:6; 7:12; 13:18};
    elseif numel(cur_elec_contact_ind) >= 15
        cnum_runs = {1:5; 6:11; 12:numel(cur_elec_contact_ind)};
    else
        cnum_runs = {1:floor(numel(cur_elec_contact_ind)/2); floor(numel(cur_elec_contact_ind)/2)+1:numel(cur_elec_contact_ind)};
    end

    for crun = 1:numel(cnum_runs)
        cnums = cnum_runs{crun};
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
        colors = colors(1:floor(height(colors) / numel(cnums)):height(colors),:);

        fg2 = figure;
        set( fg2, 'position', [ 88  100  1000  800 ] )
        tiledlayout(2,1,'TileSpacing','compact');

        global_ymin = inf;
        global_ymax = -inf;
        ax_handles = [];

        num_events = find_num_trials(subject_ID, reference, alignments);
        for alignment = alignments
            align_name = strrep(alignment, '_', ' '); 

            ax = nexttile;
            ax_handles = [ax_handles; ax];

            for ct = cnums
                cur_ct_ITPC = ITPC(subject_ID, reference, alignment, cur_letter, ct, num_events);
                plot((1:numel(cur_ct_ITPC))/Fs,cur_ct_ITPC,'linewidth',2,'color',colors(ct - cnums(1) + 1,:))
                hold on
            end
            xlim([0, numel(cur_ct_ITPC)/Fs]);
            xlabel('time from alignment (s)'); 
            ylabel('ITPC');
            title(sprintf("%s ITPC Time Series",align_name))
            legend(cur_elec_contact_names(cnums))
            %fname = sprintf("%s/%s_%s-%s.jpg", plot_out_dir, alignment, cname_st, cname_end);
            cur_ylim = ylim(ax);
            global_ymin = min(global_ymin, cur_ylim(1));
            global_ymax = max(global_ymax, cur_ylim(2));

        end % end loop through alignments
        for ax = ax_handles'
            ylim(ax, [global_ymin, global_ymax]);
        end
        cname_st = cur_elec_contact_names(cnums(end));
        cname_end = cur_elec_contact_names(cnums(1));
        fname = sprintf("%s/%s_vs_%s_%s-%s_trialct_adj.jpg", plot_out_dir, alignments(1), alignments(2), cname_st, cname_end); %all4inspections
        print(fg2,'-djpeg',fname)
        close;
    end
end % end loop through elec letters

%% Function definition

function ITPC_avg = ITPC(subject_ID, reference, alignment, cur_letter, ct, num_events)

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

ITPC_ts = zeros(max(win), 1);
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

    contact = cur_elec_contact_ind(ct);
    if strcmp(reference,'Ground')
        load(sprintf('%s/%s/separate_channel_files/%s_MAD_SES%d_ch%03d.mat',data_base_dir,subject_ID,subject_ID,sesnum,contact),"data"); 
    elseif strcmp(reference,'neighbor_average')            
        load(sprintf('%s/%s/separate_channel_files/%s/%s_MAD_SES%d_ch%03d.mat',data_base_dir,subject_ID,reference,subject_ID,sesnum,contact),"data");
    end
        
    [b,a] = butter( 4, [5 40] ./ (Fs/2) ); filtered = filtfilt( b, a, data ); data = filtered;

    rng(sesnum) % seed = sesnum
    rand_events = randperm(numel(align_times), num_events(sesnum));
    for event = rand_events
        data_padded = data(round(align_times(event)-win(1)) - Fs:round(align_times(event)+win(2)-1) + Fs);
        xgp = generalized_phase_vector(data_padded, Fs, 5); xgp = xgp(1+Fs:1+Fs+win(2)-1);
        angle_ts = angle(xgp); % time series of angles

        ITPC_ts = ITPC_ts + exp(angle_ts * 1i);
    end
end
ITPC_avg = abs(ITPC_ts / sum(num_events));

end


function num_trials = find_num_trials(subject_ID, reference, alignments)
    % Returns column vector where each element is number of trials to use for a session

    addpath('/media/Data/Human_Intracranial_MAD/analysis/TravWaves/Code')
    [~, data_base_dir, ~, ~, num_sessions, ~, ~] = tw_setup(subject_ID, reference);
    num_trials_all_align = zeros(num_sessions, numel(alignments));

    for align_num = 1:numel(alignments)
        alignment = alignments(align_num);
        inspect_alignments = ["first_unique_attribute", "second_unique_attribute", "third_unique_attribute", "fourth_unique_attribute", "inspection", "single_opt_first_inspection", "full_single_opt_info", "amount", "probability"];
        if ismember(alignment, inspect_alignments)
            align_type = "inspection";
        end
        for sesnum = 1:num_sessions   
            if strcmp(reference,'Ground')
                load(sprintf('%s/%s/%s_MAD_SES%d_Setup.mat',data_base_dir,subject_ID,subject_ID,sesnum), "filters", "trial_times", "trial_words");
            else
                load(sprintf('%s/%s/%s_MAD_SES%d_Setup_%s.mat',data_base_dir,subject_ID,subject_ID,sesnum,reference), "filters", "trial_times", "trial_words");
            end

            [align_times,~] = get_align_times(filters, trial_times, trial_words, alignment);
            align_times(isnan(align_times)) = [];

            if strcmp(subject_ID,'EMU001') && strcmp(align_type, 'inspection')
                [align_times_inspect,~] = get_align_times(filters, trial_times, trial_words, '2opt_2att_inspection');
                align_times_inspect(isnan(align_times_inspect)) = [];
                align_times = intersect(align_times, align_times_inspect);
            end
            
            num_trials_all_align(sesnum, align_num) = numel(align_times);

        end
    end

    num_trials = min(num_trials_all_align, [], 2);
end