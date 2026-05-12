function [time_jittered] = time_jitter(subject_ID, reference, iterations, contact, alignment)

% Input
% - iterations: number of time-jittered trial-avg signals to produce
% - contact: channel number from cur_elec_contact_ind corresponding to the channel to be tested

% Output
% - if iterations = 1, output the time-jittered signal
% - if iterations >> 1, output the distribution (array) of peak and trough
%   amplitudes from the generated time jittered signals

addpath('/media/Data/Human_Intracranial_MAD/analysis/TravWaves/Code')
[~, data_base_dir, ~, ~, num_sessions, Fs, ~] = tw_setup(subject_ID, reference);
if strcmp(alignment,"anticipation") || strcmp(alignment,"trialstart")
    win = [Fs/2 0];
else
    win = [0 Fs/2];
end

peak_trough_distr = [];
for iter = 1:iterations

    sum_events = zeros(max(win),1);
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
            load(sprintf('%s/%s/separate_channel_files/%s_MAD_SES%d_ch%03d.mat',data_base_dir,subject_ID,subject_ID,sesnum,contact),"data"); 
        elseif strcmp(reference,'neighbor_average')
            load(sprintf('%s/%s/separate_channel_files/%s/%s_MAD_SES%d_ch%03d.mat',data_base_dir,subject_ID,reference,subject_ID,sesnum,contact),"data");
        end

        % RNG for random time offsets
        rng(iter+sesnum, 'twister') % seed = iter+sesnum here
        rand_nums = rand(numel(align_times),1);

        % broadband filter
        [b,a] = butter( 4, [5 40] ./ (Fs/2) ); filtered = filtfilt( b, a, data ); 
        data = filtered;

        for ii = 1:numel(align_times)

            offset_sec = -0.25 + 0.5 * rand_nums(ii);
            offset_samples = offset_sec * Fs;

            data_event = data((round(align_times(ii)-win(1)+offset_samples):round(align_times(ii)+win(2)-1+offset_samples)));
            data_event = (data_event - mean(data_event))./ std(data_event); 
            data_event_padded = data(align_times(ii)-win(1)-Fs:align_times(ii)+win(2)-1+Fs);

            sum_events = sum_events + data_event;    
            sum_events_padded = sum_events_padded + data_event_padded;

        end % end loop through events
        num_events = num_events + numel(align_times);
    end % end loop through sesnums
    avg_events_jittered = sum_events ./ num_events;
    avg_events_jittered_padded = sum_events_padded ./ num_events;
    

    % % FIND PEAKS AND TROUGHS
    lp = 5; xgp_avg_jit = generalized_phase_vector( avg_events_jittered_padded, Fs, lp );
    xgp_avg_jit = xgp_avg_jit(1+Fs:1+Fs+win(2)-1);
    angle_avg_jit = angle(xgp_avg_jit);
    angle_posneg = angle_avg_jit;
    angle_posneg(angle_posneg >= 0) = 1;
    angle_posneg(angle_posneg < 0) = -1;

    peaks_jit = avg_events_jittered(diff(angle_posneg) > 0);
    troughs_jit = avg_events_jittered(diff(angle_posneg) < 0);

    % add peak and tough amplitudes to respective distributions
    peak_trough_distr = [peak_trough_distr; peaks_jit; troughs_jit];

    if iterations == 1
        time_jittered = avg_events_jittered;
    else
        time_jittered = peak_trough_distr;
    end

end % end loop through iterations

end