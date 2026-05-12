% plot_save_trial_avg.m

function [xgp_avg_allct1, avg_events_allct1, xgp_avg_allct2, avg_events_allct2, xgp_avg_allct3, avg_events_allct3, xgp_avg_allct4, avg_events_allct4] = plot_save_trial_avg(subject_ID, reference, align_type, alignments, align_names, cur_letter, time_jitter, to_plot, to_save)

% Plots to check for event-related potentials, saves trial-averaged data
% and generalized phase vector (xgp) of trial-averaged data 

% INPUTS
% - subject_ID: char array
% - reference: string or char array (only set up to handle Ground and neighbor_average)
% - align_type: character array (broad category of event: e.g., 'inspection' or 'outcome')
% - alignments, align_names: array of strings (maximum this function can handle is FOUR)
    % e.g., 
    % alignments = ["single_opt_first_inspection","full_single_opt_info"]; 
    % align_names = ["Single Option First Inspection","Full Single Option Info"];
    %
    % alignments = ["first_unique_attribute", "second_unique_attribute","third_unique_attribute","fourth_unique_attribute"]; 
    % align_names = ["Option 1 First Inspection", "Option 1 Second Inspection","Option 2 First Inspection", "Option 2 Second Inspection"];
% - cur_letter: string ("Z'", "M", "N'", etc)
% - time_jitter: boolean, true if want to plot time-jittered trial-avg signal
% - to_plot: boolean, true if want to plot and save
% - to_save: boolean, true if want to save dataa

% Outputs
% - xgp_avg_allct: Nt x Nc array where Nt is number of time points in the alignment window and Nc is the number of contacts on the specified probe
% - avg_events_allct: trial-averaged real-valued signal for all contacts on probe
% ^ one set of outputs for each alignment in alignments

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

xgp_avg_allct1 = []; avg_events_allct1 = [];
xgp_avg_allct2 = []; avg_events_allct2 = [];
xgp_avg_allct3 = []; avg_events_allct3 = [];
xgp_avg_allct4 = []; avg_events_allct4 = [];

if numel(alignments) > 4
    error('Cannot handle more than 4 alignments. Check that input is string (not char) array');
end

[~, data_base_dir, tw_out_dir, ~, num_sessions, Fs, ~] = tw_setup(subject_ID, reference);
win = [0 Fs/2];

% if strcmp(reference,'Ground')
%     load(sprintf('%s/%s/%s_MAD_SES%d_Setup.mat',data_base_dir,subject_ID,subject_ID,1));
%     %channel_ind = (1:numel(elec_name))';
% else
%     load(sprintf('%s/%s/%s_MAD_SES%d_Setup_%s.mat',data_base_dir,subject_ID,subject_ID,1,reference));
%     %elec_name = channel_name;
% end

[cur_elec_contact_ind, cur_elec_contact_names] = get_single_probe_contacts(reference, subject_ID, cur_letter);

ch_nums = 9:15; %1:numel(cur_elec_contact_ind);
colors = colors(1:floor(height(colors)/numel(ch_nums)):height(colors), :);

if to_plot
    fg1 = figure; 
    set( fg1, 'position', [ 88  100  1450  760 ] )
    if time_jitter
        t = tiledlayout(2,numel(alignments),'TileSpacing','compact');
    else
        t = tiledlayout(2,ceil(numel(alignments)/2),'TileSpacing','compact');
    end
end

for align_num = 1:numel(alignments)
    alignment = alignments(align_num);

    plot_out_dir = sprintf('%s/Plots/%s/%s/waveform/trial_avg/%s',tw_out_dir,reference,subject_ID,align_type);
    if ~exist(plot_out_dir, 'dir'), mkdir(plot_out_dir); end

    data_out_dir = sprintf('%s/Data/%s/%s/filtered_trial_avg/%s',tw_out_dir, reference,subject_ID, alignment);
    if ~exist(data_out_dir, 'dir'), mkdir(data_out_dir); end

    % preset random offsets
    rng(10, 'twister')
    all_ses_rand_nums = cell(1,num_sessions);
    for sesnum = 1:num_sessions
        if strcmp(reference,'Ground')
            load(sprintf('%s/%s/%s_MAD_SES%d_Setup.mat',data_base_dir,subject_ID,subject_ID,sesnum), "filters", "trial_times", "trial_words");
            %channel_ind = (1:numel(elec_name))';
        else
            load(sprintf('%s/%s/%s_MAD_SES%d_Setup_%s.mat',data_base_dir,subject_ID,subject_ID,sesnum,reference), "filters", "trial_times", "trial_words");
            %elec_name = channel_name;
        end
    
        [align_times,~] = get_align_times(filters, trial_times, trial_words, alignment);
    
        rand_nums = rand(numel(align_times),1);
        all_ses_rand_nums{sesnum} = rand_nums;
    end

    xgp_avg_allct = zeros(max(win),numel(cur_elec_contact_ind)); 
    avg_events_allct = zeros(max(win),numel(cur_elec_contact_ind)); 
    
    to_plot_aligned = [];
    to_plot_jittered = [];
    for cnum = 1:numel(cur_elec_contact_ind)
        contact = cur_elec_contact_ind(cnum);
        cur_name = cur_elec_contact_names(cnum);
        
        sum_events = zeros(max(win),1); sum_events_off = zeros(max(win), 1);
        sum_events_padded = zeros(max(win) + 2*Fs, 1);
        num_events = 0;
        for sesnum = 1:num_sessions      
            if strcmp(reference,'Ground')
                load(sprintf('%s/%s/%s_MAD_SES%d_Setup.mat',data_base_dir,subject_ID,subject_ID,sesnum), "filters", "trial_times", "trial_words");
                %channel_ind = (1:numel(elec_name))';
            else
                load(sprintf('%s/%s/%s_MAD_SES%d_Setup_%s.mat',data_base_dir,subject_ID,subject_ID,sesnum,reference), "filters", "trial_times", "trial_words");
                %elec_name = channel_name;
            end

            [align_times,~] = get_align_times(filters, trial_times, trial_words, alignment);

            align_times(isnan(align_times)) = [];
            align_times = round(align_times*Fs);  

            if strcmp(subject_ID,'EMU001')
                [align_times_inspect,~] = get_align_times(filters, trial_times, trial_words, '2opt_2att_inspection');

                align_times_inspect(isnan(align_times_inspect)) = [];
                align_times_inspect = round(align_times_inspect*Fs); 

                align_times = intersect(align_times, align_times_inspect);
            end
            
            if strcmp(reference,'Ground')
                load(sprintf('%s/%s/separate_channel_files/%s_MAD_SES%d_ch%03d.mat',data_base_dir,subject_ID,subject_ID,sesnum,contact), "data"); 
            elseif strcmp(reference,'neighbor_average')
                load(sprintf('%s/%s/separate_channel_files/%s/%s_MAD_SES%d_ch%03d.mat',data_base_dir,subject_ID,reference,subject_ID,sesnum,contact), "data");
            end

            [b,a] = butter( 4, [5 40] ./ (Fs/2) ); filtered = filtfilt( b, a, data ); data = filtered;
    
            for ii = 1:numel(align_times) 

                data_event = data((round(align_times(ii)-win(1)):round(align_times(ii)+win(2)-1)));
                data_event = ( data_event - mean(data_event) ) ./ std(data_event);
                data_event_padded = data(align_times(ii)-win(1)-Fs:align_times(ii)+win(2)-1+Fs);
                    % don't need to z-score padded data since we only use 
                    % this in xgp computation, and z-scoring does not affect phase
                sum_events = sum_events + data_event;    
                sum_events_padded = sum_events_padded + data_event_padded;
    
            end % end loop through events

            % rng(10, 'twister') % seed = 10; 
            %     % setting the seed within the loop ensures the same rand nums are genearted for each contact
            %     % BUT if subj has multiple sessions, each session will have the same rand offsets (bc this is within the loop through sessions)
            % rand_nums = rand(numel(align_times),1);

            for ii = 1:numel(align_times)

                rand_num = all_ses_rand_nums{sesnum}; rand_num = rand_num(ii);
                offset_sec = -0.25 + 0.5 * rand_num; % rand_nums(ii)
                offset_samples = offset_sec * Fs;

                data_event_off = data((round(align_times(ii)-win(1)+offset_samples):round(align_times(ii)+win(2)-1+offset_samples)));
                data_event_off = (data_event_off - mean(data_event_off))./ std(data_event_off);
                sum_events_off = sum_events_off + data_event_off;     

            end

            num_events = num_events + numel(align_times);

        end % end loop through sessions

        avg_events = sum_events./ num_events;
        avg_events_off = sum_events_off./ num_events;
        avg_events_padded = sum_events_padded ./ num_events;

        to_plot_aligned = [to_plot_aligned avg_events];
        to_plot_jittered = [to_plot_jittered avg_events_off];

        % Construct complex signal
        lp = 5;
        xgp_avg = generalized_phase_vector( avg_events_padded, Fs, lp );
        xgp_avg = xgp_avg(1+Fs:1+Fs+win(2)-1);

        if to_save
            fname_data_avg = sprintf("%s/centered_%s_xgp_%s_ch%03d.mat", data_out_dir, alignment, cur_name, contact);
            save(fname_data_avg, "xgp_avg")
            
            fname_data_avg = sprintf("%s/realavg_%s_%s_ch%03d.mat", data_out_dir, alignment, cur_name, contact);
            save(fname_data_avg, "avg_events")
        end

        xgp_avg_allct(:,cnum) = xgp_avg;
        avg_events_allct(:,cnum) = avg_events;

        %fprintf('Saved: contact %d\n', cnum)
        
    end % end loop through contacts

    eval(sprintf('xgp_avg_allct%d = xgp_avg_allct;',align_num));
    eval(sprintf('avg_events_allct%d = avg_events_allct;',align_num));
    
    % PLOT
    if to_plot
        nexttile;
        for cnum = 1:numel(ch_nums) %1:width(to_plot_aligned)
            plot((1:height(to_plot_aligned))/Fs,to_plot_aligned(:,ch_nums(cnum)),'linewidth',2,'color',colors(cnum,:)) % -1024:-1
            hold on
            xlim([0, 0.3])
            %xlim([0,height(to_plot_aligned)/Fs]); %xlim([0,Fs/2]); 
            xlabel('time from alignment (s)'); 
            ylabel('zscored amplitude');
            yl1 = ylim;
            title(sprintf("%s",align_names(align_num)))
            
        end
        legend(cur_elec_contact_names(ch_nums))

        if time_jitter
            nexttile;
            for cnum = 1:numel(ch_nums)
                plot((1:height(to_plot_jittered))/Fs,to_plot_jittered(:,ch_nums(cnum)),'linewidth',2,'color',colors(cnum,:)) % -1024:-1
                hold on
                xlim([0,height(to_plot_jittered)/Fs]); 
                xlabel('time from alignment (s)'); 
                ylabel('zscored amplitude');
                ylim(yl1);
                title(sprintf("%s time-jittered",align_names(align_num)))
            end
        end
    end
    
end % end loop through alignments

if to_plot
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

    if numel(alignments) == 2
        fname = sprintf("%s/allses_%s_vs_%s_elec%s.jpg", plot_out_dir, alignments(1), alignments(2),cur_letter);
    elseif numel(alignments) == 1
        fname = sprintf("%s/allses_%s_elec%s.jpg", plot_out_dir, alignments(1),cur_letter);
    elseif numel(alignments) == 4
        fname = sprintf("%s/TEST_all4inspections_elec%s.jpg", plot_out_dir,cur_letter);
    end

    print(fg1,'-djpeg',fname)
    close;
end

end 