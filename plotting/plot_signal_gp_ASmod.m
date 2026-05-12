% plot_signal_gp.m

subject_ID = 'EMU036'; % done emu036, 039, 024, 025
reference = 'Ground';
sesnum = 1; Fs = 2048;

datapath = '/media/Data/Human_Intracranial_MAD/';
data_base_dir = sprintf('%s1_formatted',datapath);
final_out_dir = sprintf('/media/Data/Human_Intracranial_MAD/analysis/Cross_Freq_Coupling/%s/%s/synchrony_data/generalized_phase',subject_ID, reference);
out_dir = sprintf('/media/Data/Human_Intracranial_MAD/analysis/Cross_Freq_Coupling/%s/%s/synchrony_data/generalized_phase',subject_ID, reference);
if ~exist(out_dir, 'dir'), mkdir(out_dir); end

load(sprintf('%s/%s/%s_MAD_SES%d_Setup.mat',data_base_dir,subject_ID,subject_ID,sesnum));

if strcmp(subject_ID, 'EMU025'), elec_name = elec_name(1:206); end

cTable = readtable(sprintf('/media/Data/Human_Intracranial_MAD/analysis/Location_Consistency/3D Coordinates/%s.tsv',subject_ID), "FileType","delimitedtext",'Delimiter', '\t');

channel_names_bs = cTable.Channel;
coord3d_all = cTable.SCS;

channel_names_list = [];
for i = 1:numel(channel_names_bs)
    cur_name = channel_names_bs{i};
    channel_names_list = [channel_names_list; convertCharsToStrings(cur_name)];
end

% Get 3d coordinates of all contacts
all_contact_coord = [];
for i = 1:numel(elec_name)
    idx = find(channel_names_list == convertCharsToStrings(elec_name{i}));
    if isempty(idx), continue; end
    coord3d = coord3d_all{idx};

    coord3d(1) = ''; coord3d(numel(coord3d)) = '';
    commas = find(coord3d == ',');
    x_coord = coord3d(1:commas(1)-1); x_coord = str2double(convertCharsToStrings(x_coord));
    y_coord = coord3d(commas(1)+1:commas(2)-1); y_coord = str2double(convertCharsToStrings(y_coord));
    z_coord = coord3d(commas(2)+1:numel(coord3d)); z_coord = str2double(convertCharsToStrings(z_coord));
    all_contact_coord = [all_contact_coord; [x_coord y_coord z_coord]];

end

% Get list of individual electrode letters
elec_letters = [];
for i = 1:numel(channel_names_bs)
    cur_name = channel_names_bs{i};
    if ~ismember(cur_name(1),elec_letters)
        elec_letters = [elec_letters; cur_name(1)];
    end
end

i = 3; % electrode iteration, can range from 1:numel(elec_letters)
cur_letter = elec_letters(i);

j = 1;
while ~contains(elec_name{j},cur_letter), j = j+1; end
cur_name = elec_name{j};
cur_elec_contact_names = []; % stores name of contacts on current electrode
cur_elec_contact_ind = []; % stores indices of contacts on curr elec
while strcmp(cur_name(1),cur_letter)
    cur_elec_contact_names = [cur_elec_contact_names; convertCharsToStrings(cur_name)];
    cur_elec_contact_ind = [cur_elec_contact_ind; j];
    
    j = j+1; cur_name = elec_name{j};
end

alignment = 'inspection';
[align_times,trial_numbers] = get_align_times(filters, trial_times, trial_words, alignment);
remove_ind = isnan(align_times);
align_times(isnan(align_times)) = [];
trial_numbers(remove_ind) = [];
align_times = round(align_times*Fs);
win = [0 1024];

event = 5; % can range from 1:numel(align_times)


% plotting
fg1 = figure; hold on; 
ax1 = gca; set( fg1, 'position', [ 88  1593  1250  420 ] )
for contact_iter = 1:10

    load(sprintf('%s/%s/separate_channel_files/%s_MAD_SES%d_ch%03d.mat', ...
        data_base_dir,subject_ID,subject_ID,sesnum,cur_elec_contact_ind(contact_iter))); 
    %data = data(round(align_times(event)-win(1)):round(align_times(event)+win(2)-1));
    
    filter_order = 4; Fs = 2048; lp = 5;
    dt = 1 / Fs; T = length(data) / Fs; time = dt:dt:T;
    
    [b,a] = butter( filter_order, [5 40] ./ (Fs/2) ); xf = filtfilt( b, a, data ); 
    
    [b,a] = butter( 4, [5 200] ./ (Fs/2) ); xw1 = filtfilt( b, a, data ); 
    [b,a] = butter( 4, [58 62] ./ (Fs/2), 'stop' ); xw1 = filtfilt( b, a, xw1 ); 
    [b,a] = butter( 4, [115 125] ./ (Fs/2), 'stop' ); xw1 = filtfilt( b, a, xw1 );
    
    colors = [
        0.80 0.90 1.00;
        0.71 0.81 0.96;
        0.63 0.73 0.92;
        0.54 0.64 0.87;
        0.45 0.56 0.83;
        0.36 0.47 0.78;
        0.27 0.39 0.74;
        0.18 0.30 0.69;
        0.09 0.22 0.65;
        0.00 0.13 0.60];
    
    xf = xf(round(align_times(event)-win(1)):round(align_times(event)+win(2)-1));
    time = time(round(align_times(event)-win(1)):round(align_times(event)+win(2)-1));

    plot(time, xf, 'linewidth', 3, 'color', colors(contact_iter,:)); 
end