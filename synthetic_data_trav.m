% Experimetation with synthetic data to produce traveling waves

%% 1. Sawtooth waves (mimic heavier right tail on single trial waveforms)

duration_samples = 1000; fs = 1000; n_samples = duration_samples + 2*fs; n_trials = 100;

decay_time1 = 0.2; decay_time2 = 0.25; decay_time3 = 0.3;
% for peak_time, generate random number between 0.4-0.6 (ch1), 0.35-0.65 (ch2), 0.3-0.7 (ch3)
ch1_jitter = [0.4, 0.6]; ch2_jitter = [0.4, 0.6]; ch3_jitter = [0.4, 0.6];
ch1_amp = 100; ch2_amp = 100; ch3_amp = 100;

ch1_data = zeros(n_trials, n_samples); % highest amplitude % each row is time series for one trial (one sawtooth impulse)
ch2_data = zeros(n_trials, n_samples);
ch3_data = zeros(n_trials, n_samples); % lowest amplitude

rand1 = ch1_jitter(1) + (ch1_jitter(2) - ch1_jitter(1) ) * rand(n_trials, 1);
rand2 = ch2_jitter(1) + (ch2_jitter(2) - ch2_jitter(1) ) * rand(n_trials, 1);
rand3 = ch3_jitter(1) + (ch3_jitter(2) - ch3_jitter(1) ) * rand(n_trials, 1);
for trial = 1:n_trials
    peak1 = rand1(trial); peak2 = rand2(trial); peak3 = rand3(trial);
    sd1 = create_sawtooth_impulse(n_samples, peak1 + 1, ch1_amp, decay_time1, fs);
    sd2 = create_sawtooth_impulse(n_samples, peak2 + 1, ch2_amp, decay_time2, fs);
    sd3 = create_sawtooth_impulse(n_samples, peak3 + 1, ch3_amp, decay_time3, fs);
    ch1_data(trial, :) = sd1;
    ch2_data(trial, :) = sd2;
    ch3_data(trial, :) = sd3;
end
avg_sd1 = mean(ch1_data, 1); avg_sd2 = mean(ch2_data, 1);  avg_sd3 = mean(ch3_data, 1); 
win = 50; % moving avg window size (samples)
avg_smooth1 = movmean(avg_sd1, win); avg_smooth2 = movmean(avg_sd2, win); avg_smooth3 = movmean(avg_sd3, win);
avg_smooth1 = avg_smooth1(1+fs:1+fs+duration_samples); avg_smooth2 = avg_smooth2(1+fs:1+fs+duration_samples); avg_smooth3 = avg_smooth3(1+fs:1+fs+duration_samples);
%[b,a] = butter( 4, [5 40] ./ (fs/2) ); 
%avg_sd1_filt = filtfilt( b, a, avg_sd1 ); avg_sd2_filt = filtfilt( b, a, avg_sd2 ); avg_sd3_filt = filtfilt( b, a, avg_sd3 ); 
%avg_sd1_filt = avg_sd1_filt(1+fs:1+fs+duration_samples); avg_sd2_filt = avg_sd2_filt(1+fs:1+fs+duration_samples); avg_sd3_filt = avg_sd3_filt(1+fs:1+fs+duration_samples);
avg_sd1_raw = avg_sd1(1+fs:1+fs+duration_samples);
avg_sd2_raw = avg_sd2(1+fs:1+fs+duration_samples);
avg_sd3_raw = avg_sd3(1+fs:1+fs+duration_samples);

fg1 = figure; % plot raw and smoothed side by side
set( fg1, 'position', [ 88  100  1450  760 ] )
tiledlayout(1,2,'TileSpacing','compact');

nexttile;
plot(avg_sd1_raw, 'LineWidth',2)
hold on
plot(avg_sd2_raw, 'LineWidth',2)
plot(avg_sd3_raw, 'LineWidth',2)
title("Full time res")
legend(["ch1", "ch2", "ch3"])

nexttile;
plot(avg_smooth1, 'LineWidth',2)
hold on
plot(avg_smooth2, 'LineWidth',2)
plot(avg_smooth3, 'LineWidth',2)
title("Smoothed")
legend(["ch1", "ch2", "ch3"])

function waveform = create_sawtooth_impulse(n_samples, peak_time, amplitude, decay_time, fs)
    % CREATE_SAWTOOTH_IMPULSE Creates a single sawtooth impulse
    %
    % Inputs:
    %   n_samples   - Total number of samples in the waveform
    %   peak_time   - Time of the peak in seconds
    %   amplitude   - Peak amplitude of the impulse
    %   decay_time  - Time constant for linear decay (time to return to zero) in seconds
    %   fs          - Sampling rate in Hz
    % Output:
    %   waveform    - Row vector of length n_samples containing the sawtooth impulse
    %
    % Example:
    %   wave = create_sawtooth_impulse(1000, 0.1, 5, 0.05, 1000);
    %   % Creates 1000 samples, peak at 100ms, amplitude 5, 50ms decay, 1kHz sampling
    
    waveform = zeros(1, n_samples);
    
    peak_sample = round(peak_time * fs);
    decay_samples = round(decay_time * fs);
    
    if peak_sample < 1 || peak_sample > n_samples
        error('Peak time is outside the waveform duration');
    end
    
    % Create sharp rise (can adjust rise_samples for steepness)
    rise_samples = max(1, round(0.01 * fs));  % 10ms rise time by default
    start_sample = max(1, peak_sample - rise_samples);
    
    % Sharp increase to peak
    rise_indices = start_sample:peak_sample;
    waveform(rise_indices) = linspace(0, amplitude, length(rise_indices));
    
    % Linear decay from peak
    end_sample = min(n_samples, peak_sample + decay_samples);
    decay_indices = peak_sample:end_sample;
    waveform(decay_indices) = linspace(amplitude, 0, length(decay_indices));
    
end