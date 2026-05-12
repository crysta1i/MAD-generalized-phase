function data_allct_alignments = plot_single_trial(subject_ID, reference, cur_letter, ch_nums, alignments, sesnum, trials)

addpath('/media/Data/Human_Intracranial_MAD/analysis/TravWaves/Code')
    
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

colors = colors(1:floor(height(colors)/numel(ch_nums))-1:height(colors), :);

[~, data_base_dir, tw_out_dir, ~, ~, Fs, ~] = tw_setup(subject_ID, reference);
save_folder = sprintf("%s/Plots/%s/%s/waveform/single_trial/Elec %s/ses%d", tw_out_dir, reference, subject_ID, cur_letter, sesnum);
if ~exist(save_folder, "dir"), mkdir(save_folder); end
win = [0 Fs/2];

[cur_elec_contact_ind, cur_elec_contact_names] = get_single_probe_contacts(reference, subject_ID, cur_letter);

if strcmp(reference,'Ground')
    load(sprintf('%s/%s/%s_MAD_SES%d_Setup.mat',data_base_dir,subject_ID,subject_ID,sesnum), "filters", "trial_times", "trial_words");
    %channel_ind = (1:numel(elec_name))';
else
    load(sprintf('%s/%s/%s_MAD_SES%d_Setup_%s.mat',data_base_dir,subject_ID,subject_ID,sesnum,reference), "filters", "trial_times", "trial_words");
    %elec_name = channel_name;
end

% reconcile number of trials across alignments
min_num_trials = nan;
for alignment = alignments
    [align_times,~] = get_align_times(filters, trial_times, trial_words, alignment);
    min_num_trials = min_num_trials + numel(align_times);
end
trials(trials > min_num_trials) = [];

all_align_times = cell(1, numel(alignments));
for align_num = 1:numel(alignments)
    alignment = alignments(align_num);
    [align_times,~] = get_align_times(filters, trial_times, trial_words, alignment);
    align_times(isnan(align_times)) = [];
    align_times = round(align_times*Fs);  

    if strcmp(subject_ID,'EMU001')
        [align_times_inspect,~] = get_align_times(filters, trial_times, trial_words, '2opt_2att_inspection');

        align_times_inspect(isnan(align_times_inspect)) = [];
        align_times_inspect = round(align_times_inspect*Fs); 

        align_times = intersect(align_times, align_times_inspect);
    end

    all_align_times{align_num} = align_times;
end

data_allct_alignments = cell(numel(ch_nums),1); % each entry of cell is a 3D array for one contact: trial x timepts x alignments
for ct = 1:numel(ch_nums)
    contact = cur_elec_contact_ind(ch_nums(ct));
    if strcmp(reference,'Ground')
        load(sprintf('%s/%s/separate_channel_files/%s_MAD_SES%d_ch%03d.mat',data_base_dir,subject_ID,subject_ID,sesnum,contact), "data"); 
    elseif strcmp(reference,'neighbor_average')
        load(sprintf('%s/%s/separate_channel_files/%s/%s_MAD_SES%d_ch%03d.mat',data_base_dir,subject_ID,reference,subject_ID,sesnum,contact), "data");
    end

    cur_ct_data = zeros(numel(trials), abs(win(2)-win(1)), numel(alignments));
    [b,a] = butter( 4, [5 40] ./ (Fs/2) ); filtered = filtfilt( b, a, data ); data = filtered;

    for ii = 1:numel(trials)
        cur_trial_data = [];
        event = trials(ii);
        for align_num = 1:numel(alignments)
            align_times = all_align_times{align_num};
            data_event = data((round(align_times(event)-win(1)):round(align_times(event)+win(2)-1)));
            data_event = ( data_event - mean(data_event) ) ./ std(data_event);
            %data_event_padded = data(align_times(ii)-win(1)-Fs:align_times(ii)+win(2)-1+Fs);
            cur_trial_data = [cur_trial_data data_event]; % timepts x alignment 
        end
        cur_ct_data(ii,:,:) = cur_trial_data;
    end
    data_allct_alignments{ct} = cur_ct_data;
end % end loop thru contacts

for event = 1:numel(trials)
    fg1 = figure; 
    set( fg1, 'position', [ 88  100  1450  760 ] )
    t = tiledlayout(2,ceil(numel(alignments)/2),'TileSpacing','compact');

    for align_num = 1:numel(alignments)
        alignment = alignments(align_num);
        align_name = strrep(alignment, '_', ' ');
        nexttile;
        ct_names = [];
        for cnum = 1:numel(ch_nums)
            ct_names = [ct_names; cur_elec_contact_names(ch_nums(cnum))];
            cur_ct_data = data_allct_alignments{cnum};
            plot((1:size(cur_ct_data,2))/Fs,cur_ct_data(event, :, align_num),'linewidth',2,'color',colors(cnum,:)) % -1024:-1
            hold on
        end
        xlim([0 0.3]) %xlim([0,size(cur_ct_data,2)/Fs]); %xlim([0,Fs/2]); 
        xlabel('time from alignment (s)'); 
        ylabel('zscored amplitude');
        % if align_num == 1
        %     yl1 = ylim;
        % else
        %     ylim(yl1);
        % end
        title(sprintf("%s",align_name))
        legend(ct_names)
    end

    axObjects = findobj(t, 'Type', 'axes');

    allYLim = zeros(length(axObjects), 2);
    allRanges = zeros(length(axObjects), 1);
    for i = 1:length(axObjects)
        allYLim(i, :) = axObjects(i).YLim;               % Get current [ymin, ymax]
        allRanges(i) = allYLim(i, 2) - allYLim(i, 1);    % Calculate range span
    end
    [~, maxIdx] = max(allRanges);
    targetLimits = allYLim(maxIdx, :);
    set(axObjects, 'YLim', targetLimits);

    % SAVE
    name1 = cur_elec_contact_names(ch_nums(end)); name2 = cur_elec_contact_names(ch_nums(1));
    if numel(alignments) == 2
        fname = sprintf("%s/allses_%s_vs_%s_%s-%s_trial%d.jpg", save_folder, alignments(1), alignments(2), name1, name2, trials(event));
    elseif numel(alignments) == 1
        fname = sprintf("%s/allses_%s_%s-%s_trial%d.jpg", save_folder, alignments(1), name1, name2, trials(event));
    elseif numel(alignments) == 4
        fname = sprintf("%s/all4inspections_%s-%s_trial%d.jpg", save_folder, name1, name2, trials(event));
    end

    print(fg1,'-djpeg',fname)
    close;

end % end loop through trials to plot

end