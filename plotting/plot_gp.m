%% plot (2D and 3D) coordinates of (outer) contacts, colored by phase

%% 1.1 Get (outer) contact coordinates

subject_ID = 'EMU036'; % done emu036, 039 emu024, 025
reference = 'Ground';
sesnum = 1; Fs = 2048;

datapath = '/media/Data/Human_Intracranial_MAD/';
data_base_dir = sprintf('%s1_formatted',datapath);
final_out_dir = sprintf('/media/Data/Human_Intracranial_MAD/analysis/Cross_Freq_Coupling/%s/%s/synchrony_data/generalized_phase',subject_ID, reference);
out_dir = sprintf('/media/Data/Human_Intracranial_MAD/analysis/Cross_Freq_Coupling/%s/%s/synchrony_data/generalized_phase',subject_ID, reference);
if ~exist(out_dir, 'dir'), mkdir(out_dir); end

load(sprintf('%s/%s/%s_MAD_SES%d_Setup.mat',data_base_dir,subject_ID,subject_ID,sesnum));

%if strcmp(subject_ID, 'EMU025'), elec_name = elec_name(1:206); end

% load brainstorm file and outer contacts
cTable = readtable(sprintf('/media/Data/Human_Intracranial_MAD/analysis/Location_Consistency/3D Coordinates/%s.tsv',subject_ID), "FileType","delimitedtext",'Delimiter', '\t');
load(sprintf('%s/extracted_phase/outer_contacts/outer_contact_names.mat',out_dir)); % outer_contacts
load(sprintf('%s/extracted_phase/outer_contacts/outer_contacts_ind.mat',out_dir));

% if strcmp(subject_ID, 'EMU025') 
%     outer_contacts(5) = []; % remove M15 (right hemisphere electrode)
%     outer_contacts_ind(5) = [];
% end

channel_names_bs = cTable.Channel; % extract "Channels" variable in table
coord3d_all = cTable.SCS;

channel_names_list = [];
for i = 1:numel(channel_names_bs)
    cur_name = channel_names_bs{i};
    channel_names_list = [channel_names_list; convertCharsToStrings(cur_name)];
end

if exist(sprintf('%s/extracted_phase/outer_contacts/outer_contact_coord.mat',out_dir),'file')
    load(sprintf('%s/extracted_phase/outer_contacts/outer_contact_coord.mat',out_dir));
    % if strcmp(subject_ID, 'EMU025') 
    %     outer_contact_coord(5,:) = []; % remove M15 (right hemisphere electrode)
    % end
else
    outer_contact_coord = []; % store coordinates of all outer contacts

    % NOTE: contact coordinates are stored in the order of contacts in elec_name
    for i = 1:numel(outer_contacts)
        outer_idx = find(channel_names_list == outer_contacts(i)); %outer contacts values are strings
        coord3d = coord3d_all(outer_idx);
        coord3d = coord3d{1};
        coord3d(1) = ''; coord3d(numel(coord3d)) = '';
        commas = find(coord3d == ',');
        x_coord = coord3d(1:commas(1)-1); x_coord = str2double(convertCharsToStrings(x_coord));
        y_coord = coord3d(commas(1)+1:commas(2)-1); y_coord = str2double(convertCharsToStrings(y_coord));
        z_coord = coord3d(commas(2)+1:numel(coord3d)); z_coord = str2double(convertCharsToStrings(z_coord));
        outer_contact_coord = [outer_contact_coord; [x_coord y_coord z_coord]]; %{outer_contact_coord; channel_names{outer_idx} t{outer_idx, "SCS"}};
    
    end

    % for EMU036: y-coordinate has smallest variance and smallest range -> drop y-coord
    % (plot only 1st and 3rd columns of outer_contact_coord)

    % save all three columns of coordinates
    fname_coord = sprintf('%s/extracted_phase/outer_contacts/outer_contact_coord.mat',out_dir);
    save(fname_coord, "outer_contact_coord")
end

%% 1.2 Get all contact coordinates
% Format so that each column is a different dimension (x-coord, y-coord, z-coord)

% NOTE: contact coordinates are store3.5d in the order of contacts in elec_name
all_contact_coord = [];
for i = 1:numel(elec_name)
    idx = find(channel_names_list == convertCharsToStrings(elec_name{i})); %outer contacts values are strings
    if isempty(idx), continue; end
    coord3d = coord3d_all{idx};

    coord3d(1) = ''; coord3d(numel(coord3d)) = '';
    commas = find(coord3d == ',');
    x_coord = coord3d(1:commas(1)-1); x_coord = str2double(convertCharsToStrings(x_coord));
    y_coord = coord3d(commas(1)+1:commas(2)-1); y_coord = str2double(convertCharsToStrings(y_coord));
    z_coord = coord3d(commas(2)+1:numel(coord3d)); z_coord = str2double(convertCharsToStrings(z_coord));
    all_contact_coord = [all_contact_coord; [x_coord y_coord z_coord]];

end

%% 2. PLOT: static individual plots
addpath("/media/Data/Human_Intracranial_MAD/_toolbox/spectral-analysis")

% 2D - outer contacts only (project 3d coordinates onto 2d plane)
% load table of coordinates, outer_contacts, and gp (data) for each subject
reference = 'Ground'; 
subject_IDs = {'EMU001','EMU024', 'EMU025','EMU030', 'EMU036', 'EMU037','EMU038','EMU039','EMU041','EMU051'};
for subject_num = 3 %[2 3 4 5 6 7 8 10] 
    %
    subject_ID = subject_IDs{subject_num};
    % TODO: create function to replace switch-case block (in all functions)
    switch subject_ID
        case 'EMU001'
            num_sessions = 3;
            Fs = 1000;
        case 'EMU024'
            num_sessions = 3;
            Fs = 2048;
        case 'EMU025'
            num_sessions = 2;
            Fs = 2048;
        case 'EMU030'
            num_sessions = 2;
            Fs = 2048;
        case 'EMU036'
            num_sessions = 4;
            Fs = 2048;
        case 'EMU037'
            num_sessions = 4;
            Fs = 2048;
        case 'EMU038'
            num_sessions = 1;
            Fs = 2048;
        case 'EMU039'
            num_sessions = 4;
            Fs = 2048;
        case 'EMU041'
            num_sessions = 9;
            Fs = 2048;
        case 'EMU051'
            num_sessions = 1;
            Fs = 2048;
    end
    
    % TEMP vars
    alignment = 'inspection';
    sesnum = 1;

    %window = 'feedback';           
    % set length for aligned data window
    time_window_seconds = 1;
    % set sampling rate-dependent spectral parameters
    time_window_length = time_window_seconds*Fs;
    win = [0 time_window_length];

    %paths
    datapath = '/media/Data/Human_Intracranial_MAD/';
    data_base_dir = sprintf('%s1_formatted',datapath);
    final_out_dir = sprintf('/media/Data/Human_Intracranial_MAD/analysis/Cross_Freq_Coupling/%s/%s/synchrony_data/generalized_phase',subject_ID, reference);
    out_dir = sprintf('/media/Data/Human_Intracranial_MAD/analysis/Cross_Freq_Coupling/%s/%s/synchrony_data/generalized_phase',subject_ID, reference);
    if ~exist(out_dir, 'dir'), mkdir(out_dir); end
    
    % load outer_contacts_ind and outer_contacts
    %load(sprintf("%s/extracted_phase/outer_contacts/outer_contacts_ind.mat", out_dir));
    %load(sprintf("%s/extracted_phase/outer_contacts/outer_contact_names.mat",out_dir)); % TEMP

    % load 3d coordinates
    fname_coord = sprintf('%s/extracted_phase/outer_contacts/outer_contact_coord.mat',out_dir);
    load(fname_coord, "outer_contact_coord")
    
    % load ind_regions_all and regions_all - here, reference = neighbor_average 
    % because these variables were saved in neighbor_avg folder
    fname = sprintf('/media/Data/Human_Intracranial_MAD/analysis/Cross_Freq_Coupling/%s/%s/synchrony_data/index_regions_relevantChannels.mat',subject_ID, 'neighbor_average');
    load(fname);
    load(sprintf('/media/Data/Human_Intracranial_MAD/analysis/Cross_Freq_Coupling/%s/%s/synchrony_data/regions_relevantChannels.mat',subject_ID, 'neighbor_average'));   

    % load setup files
    load(sprintf('%s/%s/%s_MAD_SES%d_Setup.mat',data_base_dir,subject_ID,subject_ID,sesnum));

    % load all contact 3d coordinates
    %cTable = readtable("/media/Data/Human_Intracranial_MAD/analysis/Location_Consistency/3D Coordinates/EMU036.tsv", "FileType","delimitedtext",'Delimiter', '\t');
    %channel_names = cTable.Channel; % extract "Channels" variable in table
    xgp_to_plot = [];
    for time_pt_iter = 390 %TODO [3 6 9 12 15 18 21 24]

        for i = 1:numel(outer_contacts)
            % outer_contact_coord = []; % store coordinates of all outer contacts
            % outer_idx = find(convertCharsToStrings(channel_names) == outer_contacts(i)); %outer contacts values are strings
            % coord3d = cTable{outer_idx, "SCS"}{1};
            % coord3d(1) = ''; coord3d(numel(coord3d)) = '';
            % commas = find(coord3d == ',');
            % x_coord = coord3d(1:commas(1)-1);
            % y_coord = coord3d(commas(1)+1:commas(2)-1);
            % z_coord = coord3d(commas(2)+1:numel(coord3d));
            % outer_contact_coord = [outer_contact_coord; [x_coord y_coord z_coord]]; %{outer_contact_coord; channel_names{outer_idx} t{outer_idx, "SCS"}};
    
    
            % load xgp 

            ind_cur = outer_contacts_ind(i); % changed from: elec_ind(outer_contacts_ind(i));
    
            % find where ind_regions_all == outer_contacts_ind(i), and get corresponding regions_all element
            ind_idx = find(ind_regions_all == outer_contacts_ind(i));
            region_name = regions_all(ind_idx);
    
            load(sprintf('%s/extracted_phase/outer_contacts/xgp_ses%02d_ch%03d.mat',final_out_dir,sesnum,ind_cur));
            
            % extract only data from window of interest (TBD/TODO)
            [align_times,trial_numbers] = get_align_times(filters, trial_times, trial_words, alignment);
    
            % remove NaN values
            remove_ind = isnan(align_times);
            align_times(isnan(align_times)) = [];
            trial_numbers(remove_ind) = []; 
            
            % convert alignment times from seconds to samples
            align_times = round(align_times*Fs);
    
            for event = 4 %:numel(align_times) % only 1st event
                if align_times(event)+win(2)-1>size(xgp,1) || align_times(event)-win(1)-1<1
                    continue
                else
                    xgp_window = xgp((round(align_times(event)-win(1)):round(align_times(event)+win(2)-1)));
                    xgp_pt = xgp_window(time_pt_iter); % TBD/TODO: select time points that are separated by a few samples?
                    xgp_to_plot = [xgp_to_plot; xgp_pt];
                end
            end
        
        end % end loop through outer contacts

        % project onto 2D (get rid of one of the coordinates)

        % RUN THIS IN CMD WINDOW compute variance/range of each of the columns in coord3d
        % varX = var(outer_contact_coord(:,1)); varY = var(outer_contact_coord(:,2)); varZ = var(outer_contact_coord(:,3));
        % rangeX = range(outer_contact_coord(:,1)); rangeY = range(outer_contact_coord(:,2)); rangeZ = range(outer_contact_coord(:,3));
        % 
        % if min([rangeX rangeY rangeZ]) == rangeX
        %     %elimDim = 'x'; % ONLY USE y_coord and z_coord
        % elseif min([rangeX rangeY rangeZ]) == rangeY
        %     %elimDim = 'y'; % ONLY USE x_coord and z_coord
        % elseif min([rang        eX rangeY rangeZ]) == rangeZ
        %     %elimDim = 'z'; % ONLY USE x_coord and y_coord
        % end
    
        % PLOT X and Z coordinates, (circular) color bar for different phase values
        % determine color based on xgp_to_plot.

        % fname_coord = sprintf('%s/extracted_phase/outer_contacts/outer_contact_coord.mat',out_dir);
        % load(fname_coord, "outer_contact_coord")

        %fg1 = figure; hold on; ax1 = gca; set( fg1, 'position', [ 88  1593  1250  420 ] ) % TODO
        %plot(outer_contact_coord(:,1),outer_contact_coord(:,3),"Marker",'o',"MarkerSize",9) 
        %map = colorcet( 'C2' ); colormap( circshift( map, [ 28, 0 ] ) )


        phase_angles = angle(xgp_to_plot);
        %cmap = colorcet('C2'); % colorcet('C2','N',500); % 500 different colors 
        cmap = blue_cyclic_colormap(256);
        num_colors = size(cmap, 1);
        ph_min = min(phase_angles); ph_max = max(phase_angles);
        c_norm_angles = (phase_angles - ph_min) / (ph_max - ph_min); % convert from [-pi pi] to [0 1] % TODO: check that this works!!
        c_indices = round(c_norm_angles * (num_colors - 1)) + 1; % convert from [0 1] to [0 256], since 256 different colors in cmap
        
        % select colors in cmap that correspond to the value of the phase angles
        c_angles = cmap(c_indices,:);

        fg1 = figure; 
        scatter(outer_contact_coord(:,1), outer_contact_coord(:,3), 200, c_angles, 'filled') % TODO: 200 is marker size
        colormap(cmap) % Apply cyclic colormap
        colorbar
        clim([-pi pi]) % Set colorbar range to match phase angles % changed from caxis to clim
        xlim([10, 80]) % TODO different for each subj
        ylim([20, 125]) % TODO
        xlabel("x coordinate (mm)")
        ylabel("z coordinate (mm)")

        hold on
        % for i = 1:numel(outer_contacts)
        %     text(x(i) + 0.05, y(i), outer_contacts(i), 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'left');
        % end

        fname = sprintf('%s/surface_waves_fig_sfn25.jpg',final_out_dir);
        print(fg1,'-djpeg',fname)
        %close;

    end

    % TODO LATER: project onto nearest 2D plane (use all three coordinates, change of axis/basis)

end


%% ** Synthetic Data

time = 1:300;
clear angle;

% Data setup: Random noise
data_noise = zeros(numel(time),height(outer_contact_coord));
for i = 1:height(outer_contact_coord)
    rand_angles = 2 * pi * rand(1,800);
    data_noise(:,i) = rand_angles'; 
end

data = data_noise;

% Data setup: contacts each move smoothly between 0-2pi but not synchronized
data_rand = zeros(numel(time),height(outer_contact_coord)); 
rand_angles = 2 * pi * rand(1,height(outer_contact_coord)); % initialize random start angle for each contact
for i = 1:height(outer_contact_coord)
    for j = 1:numel(time)
        data_rand(j,i) = sin(rand_angles(i) + 0.3*j); 
    end
end
phase_rand = asin(data_rand);

data = phase_rand;

% Data setup: Synchronized (standing wave)
    % can imagine what this would look like


% Data setup: Wavefront moving
    % horizontal: sin { (x + (0 - min_x)) * (7/ (max_x - min_x) ) + (0.1 * t) }
synth_horiz = zeros(numel(time),height(outer_contact_coord));
for i = 1:numel(time)
    for j = 1:height(outer_contact_coord)
        x_coor = outer_contact_coord(j,1);
        norm_x = (x_coor + 45) / 90; % max_x - min_x = 90, % 0 - min_x = 45
        t = i;
        synth_horiz(i, j) = sin(2*pi*norm_x + ((3*2*pi)/300)*t);
        %synth_horiz(i, j) = sin( (x_coor + 45)*( 2*pi /(90)) + 0.3*i );
    end
end
phase_horiz = asin(synth_horiz); 
data = phase_horiz;

    % vertical: sin { (y - min_y) * (7/ (max_y - min_y) ) + (0.1 * t) } % t = i
synth_vert = zeros(numel(time),height(outer_contact_coord));
for i = 1:numel(time)
    for j = 1:height(outer_contact_coord)
        z_coor = outer_contact_coord(j,3); % 3rd column
        norm_z = (z_coor - 20) / 70;
        t = i;
        synth_vert(i, j) = sin(2*pi*norm_z + (6*pi/300)*t);
        %synth_vert(i, j) = sin((z_coor - 20)*(2*pi/(70)) + 0.3*i); % TODO
    end
end
phase_vert = asin(synth_vert);
data = phase_vert;


%% ALL DATA: create visualization (color gradient):

addpath('/media/Data/Human_Intracranial_MAD/analysis/PAC_code/CF_Coupling/generalized phase');

cmap = blue_cyclic_colormap(256); % using custom blue colormap
%cmap = colorcet('C2'); % colorcet('C2','N',500); % 500 different colors
    % 'C2' = cyclic rainbow, 'C5' = monochrome
num_colors = size(cmap, 1);

colors_ind_all = zeros(height(data),width(data)); % each column = 1 outer contact, each row = 1 timept

for i = 1:height(data) % iterate through time points
    
    cur_timept_angles = data(i,:);

    ph_min = min(cur_timept_angles); ph_max = max(cur_timept_angles);
    c_norm_angles = (cur_timept_angles - ph_min) / (ph_max - ph_min); % convert from [-pi pi] to [0 1]
    c_indices = round(c_norm_angles * (num_colors - 1)) + 1; % convert from [0 1] to [0 256] (256 colors in cmap)
    colors_ind_all(i,:) = c_indices; % use these indices to access the cmap variable

end

% Each row of colors_ind_all = one frame in animation
    % initialize first frame
figure;
cur_ind_colors = colors_ind_all(1,:); 
cur_colors = cmap(cur_ind_colors,:);
%cScatter = scatter3(outer_contact_coord(:,1), outer_contact_coord(:,2), outer_contact_coord(:,3), 200, cur_colors, 'filled'); % 200 is marker size
%cScatter = scatter(outer_contact_coord(:,1), outer_contact_coord(:,3), 200, cur_colors, 'filled');
cScatter = scatter(outer_contact_coord(:,1), 30*ones(size(outer_contact_coord(:,3))), 200, cur_colors, 'filled');
colormap(cmap) % Apply cyclic colormap
colorbar
clim([-pi pi]) % Set colorbar range to match phase angles
xlim([-45, 45]) % TODO different for each subj
ylim([20, 90])

xlabel('X') 
ylabel('Z')

% rest of the frames
for time_pt = 2:height(data)
    cur_ind_colors = colors_ind_all(time_pt,:); 
    cur_colors = cmap(cur_ind_colors,:);
    
    set(cScatter, 'CData', cur_colors);
    title(sprintf('Time Step: %d', time_pt));
    drawnow;
    pause(0.18);

end

%% ALL DATA: create binary-color visualization

wavefront = zeros(size(data));
for i = 1:height(data)
    wavefront(i,:) = abs(data(i,:)) < 0.4;
end

cur_colors = zeros(15,3);
cur_wavefront = wavefront(1,:);
for i = 1:numel(cur_wavefront)
    if cur_wavefront(i) == 1
        cur_colors(i,:) = [1 0 0];
    else
        cur_colors(i,:) = [0 0 1];
    end
end

figure;
%cScatter = scatter(outer_contact_coord(:,1), 30*ones(size(outer_contact_coord(:,3))), 200, cur_colors, 'filled');
cScatter = scatter(outer_contact_coord(:,1), outer_contact_coord(:,3), 200, cur_colors, 'filled');
xlim([-45, 45]) 
ylim([20, 90])
xlabel('X') 
ylabel('Z')

for time_pt = 2:height(data)
    cur_wavefront = wavefront(time_pt,:);
    for i = 1:numel(cur_wavefront)
        if cur_wavefront(i) == 1
            cur_colors(i,:) = [1 0 0];
        else
            cur_colors(i,:) = [0 0 1];
        end
    end
    
    set(cScatter, 'CData', cur_colors);
    title(sprintf('Time Step: %d', time_pt));
    drawnow;
    pause(0.18);

end

%% Synthetic traveling wave based on chatGPT function

x = outer_contact_coord(:,1);
y = outer_contact_coord(:,3);

% Parameters
spatial_freq = 1/40;    % cycles/unit
speed = 5;              % units/second
direction = [1, 0];     % horizontal
direction = [0, 1];     % vertical
direction = [-1, 1];     % diagonal
duration = 50;          % seconds
fps = 10;               % frames per second

% Animate
tot_wave = chatGPT_animate_traveling_wave(x, y, spatial_freq, speed, direction, duration, fps);


%% 3.0 REAL DATA - Create 2D animation of (outer contacts) scatterplots

trial_batch_analysis = 0; % TODO change this manually (= 1 if analyzing multiple inspection events, =0 if individual)

%subject_ID = 'EMU036'; sesnum = 1; reference = 'Ground'; Fs = 2048;
alignments = {'outcome','inspection','response','trialstart'};
%window_nums = {1,2}; % 1 = pre-alignment, 2 = post-alignment (e.g., pre-outcome = anticipation, post-outcome = feedback)

final_out_dir = sprintf('/media/Data/Human_Intracranial_MAD/analysis/Cross_Freq_Coupling/%s/%s/synchrony_data/generalized_phase',subject_ID, reference);
out_dir = sprintf('/media/Data/Human_Intracranial_MAD/analysis/Cross_Freq_Coupling/%s/%s/synchrony_data/generalized_phase',subject_ID, reference);
if ~exist(out_dir, 'dir'), mkdir(out_dir); end

% load outer_contact_coord, outer_contacts, and outer_contacts_ind depends only on subject_ID
fname_coord = sprintf('%s/extracted_phase/outer_contacts/outer_contact_coord.mat',out_dir);
load(fname_coord, "outer_contact_coord")
load(sprintf("%s/extracted_phase/outer_contacts/outer_contacts_ind.mat", out_dir));
load(sprintf("%s/extracted_phase/outer_contacts/outer_contact_names.mat",out_dir));

% if strcmp(subject_ID, 'EMU025') 
%     outer_contacts(5) = []; % remove M15 (right hemisphere electrode)
%     outer_contacts_ind(5) = [];
%     outer_contact_coord(5,:) = [];
% end

% load ind_regions_all and regions_all
% NOTE: here, reference = neighbor_average because these variables were saved in neighbor_avg folder
% fname = sprintf('/media/Data/Human_Intracranial_MAD/analysis/Cross_Freq_Coupling/%s/%s/synchrony_data/index_regions_relevantChannels.mat',subject_ID, 'neighbor_average');
% load(fname);
% load(sprintf('/media/Data/Human_Intracranial_MAD/analysis/Cross_Freq_Coupling/%s/%s/synchrony_data/regions_relevantChannels.mat',subject_ID, 'neighbor_average'));   

% load setup files
load(sprintf('%s/%s/%s_MAD_SES%d_Setup.mat',data_base_dir,subject_ID,subject_ID,sesnum));

xgp_full_res = [];
xgp_full_extra = [];
xgp_all_contacts = []; % each column = 1 (outer) contact, each row = 1 timept

for align_iter = 2 % CHANGE THIS MANUALLY
    alignment = alignments{align_iter};

    switch alignment
        case 'outcome'
            windows = {'anticipation','feedback'};
            window_num = 2; %TODO: manually change number
            window = windows{window_num}; % feedback
            time_window_seconds = 1; % set length for aligned data window
            time_window_length = time_window_seconds*Fs; % set sampling rate-dependent spectral parameters
            if window_num == 1
                win = [time_window_length 0];
            else
                win = [0 time_window_length];
            end
        case 'inspection'
            % only post-alignment period matters
            window = 'post-inspection';
            time_window_seconds = 0.5; % TODO: get_time_between_inspections
            time_window_length = time_window_seconds*Fs; 
            win = [0 time_window_length];

            [time_btwn_inspections, last_inspection_to_response] = filter_inspections(subject_ID, sesnum);

            inspect_types = {'win_domain','loss_domain','amount','probability','novel_attribute_inspection','repeated_inspection',...
                'first_unique_attribute','second_unique_attribute','third_unique_attribute','fourth_unique_attribute'};

            % alignment2 = subalignment (of alignment = inspection)
            alignment2 = 'novel_attribute_inspection'; % TODO: manually change this 

        case 'response'
            windows = {'pre-response','post-response'};
            window = windows{1}; % TODO: manually change number
            % TODO: time_window_length differs depending on pre- or post-alignment
            time_window_seconds = 0.5;
            time_window_length = time_window_seconds*Fs; 
            if window_num == 1
                win = [time_window_length 0];
            else
                win = [0 time_window_length];
            end
        case 'trialstart' 
            window = 'pre-trial'; % to analyze inter-trial intervals!
            time_window_seconds = 0.5;
            time_window_length = time_window_seconds*Fs;
            win = [time_window_length 0];
    end

    addpath('/media/Data/Human_Intracranial_MAD/_toolbox/spectral-analysis');
    [align_times,trial_numbers] = get_align_times(filters, trial_times, trial_words, alignment);

    remove_ind = isnan(align_times);
    align_times(isnan(align_times)) = [];
    trial_numbers(remove_ind) = [];
    
    align_times = round(align_times*Fs); % convert alignment times from seconds to samples


    for contact_iter = 1:numel(outer_contacts)

        if ~ismember(outer_contacts(contact_iter), channel_names_bs)
            continue;
        end

        ind_idx = outer_contacts_ind(contact_iter);
        
        load(sprintf('%s/extracted_phase/xgp_%s_ses%02d_ch%03d.mat',final_out_dir,elec_name{ind_idx},sesnum,ind_idx));
        % variable loaded: xgp

        if trial_batch_analysis == 1
            
            trial_num = 3; % TODO change this manually
            align_times_trial = align_times(trial_numbers == trial_num); % get all align_times (inspections) associated with this trial number
            last_idx_trial = find(trial_numbers == trial_num, 1, "last");
            time_window_samples = align_times(last_idx_trial + 1) - align_times_trial(1);
            win = [0 time_window_samples];
            xgp_window = xgp((round(align_times_trial(1)-win(1)):round(align_times_trial(1)+win(2)-1)));
            xgp_to_plot = xgp_window(3:3:numel(xgp_window));

        else
            event = 2; % TODO change this manually 
            
            % TODO: THIS DOES NOT WORK if targeting specific type of inspection
            % (only works when looking at all inspection events as a whole)
                % if (strcmp(alignment, 'response') && strcmp(window, 'pre-response'))
                %     time_window_samples = last_inspection_to_response(event);
                %     win = [time_window_samples 0]; 
                % elseif strcmp(alignment, 'inspection') && trial_numbers(event + 1) ~= trial_numbers(event)
                %     % if this event is the last inspection in a trial, use time between last inspection and response
                %     time_window_samples = last_inspection_to_response(trial_numbers(event));
                %     win = [0 time_window_samples];
                % elseif strcmp(alignment, 'inspection')
                %     time_window_samples = time_btwn_inspections(event);
                %     win = [0 time_window_samples]; % post-inspection
                % end
    
            xgp_window = xgp((round(align_times(event)-win(1)):round(align_times(event)+win(2)-1)));
            xgp_window_extra = xgp((round(align_times(event)-win(1)):round(align_times(event)+2*win(2)))); 
            xgp_to_plot = xgp_window(3:3:numel(xgp_window));
    
        end 

        xgp_full_res = [xgp_full_res xgp_window];
        xgp_full_extra = [xgp_full_extra xgp_window_extra];
        xgp_all_contacts = [xgp_all_contacts xgp_to_plot]; % HORIZONTAL
       

    end % end loop through outer contacts

end
clear angle;
angles_full_extra = angle(xgp_full_extra);
angles_full_res = angles_full_extra(1:height(xgp_full_res),:);
angles_all_contacts = angle(xgp_all_contacts);

% 1. loop through time samples to create matrix of c_indices (repeat for additional alignment points)
    % normalize each ROW of angles_all_contacts
addpath('/media/Data/Human_Intracranial_MAD/analysis/PAC_code/CF_Coupling/generalized phase');
%cmap = colorcet('C5'); % colorcet('C2','N',500); % 500 different colors
cmap = blue_cyclic_colormap(256);
% 'C2' is cyclic rainbow, 'C5' is monochrome
num_colors = size(cmap, 1);

colors_ind_all = zeros(height(angles_all_contacts),width(angles_all_contacts)); % each column = 1 outer contact, each row = 1 timept

for i = 1:height(angles_all_contacts) % iterate through time points
    
    cur_timept_angles = angles_all_contacts(i,:);

    ph_min = min(cur_timept_angles); ph_max = max(cur_timept_angles);
    c_norm_angles = (cur_timept_angles - ph_min) / (ph_max - ph_min); % convert from [-pi pi] to [0 1] % TODO: check that this works!!
    c_indices = round(c_norm_angles * (num_colors - 1)) + 1; % convert from [0 1] to [0 256], since 256 different colors in cmap
    colors_ind_all(i,:) = c_indices; % use these indices to access the cmap variable
    
    %c_angles = cmap(c_indices,:); % these are the COLORS that will be used as argument in scatter() TODO: rename this to c_colors
    %colors_all = [colors_all; c_angles]; 

end

% 2. loop through time samples again to generate colored scatter plots BASED ON matrix colors_all
% plot row-by-row of angles_all_contacts/colors_all - each row = 1 frame
% create movie here

% initialize first frame
figure('Position', [300, 50, 500, 400]);
cur_ind_colors = colors_ind_all(1,:); 
cur_colors = cmap(cur_ind_colors,:);
cScatter = scatter(outer_contact_coord(:,1), outer_contact_coord(:,3), 200, cur_colors, 'filled');
colormap(cmap) % Apply cyclic colormap
colorbar
clim([-pi pi]) % Set colorbar range to match phase angles % changed from caxis to clim

if strcmp(subject_ID, 'EMU036')
    xlim([-45, 45])
    ylim([20, 90])
elseif strcmp(subject_ID, 'EMU024')
    xlim([-30, 45])
    ylim([5, 75])
elseif strcmp(subject_ID, 'EMU025')
    xlim([10, 85])
    ylim([20, 125])
elseif strcmp(subject_ID,'EMU030')
    xlim([-25, 45])
    ylim([20, 90])
end

xlabel('X') 
ylabel('Z')

% rest of the frames
for time_pt = 2:height(angles_all_contacts)
    cur_ind_colors = colors_ind_all(time_pt,:); %cur_angles = angles_all_contacts(time_pt,:);
    cur_colors = cmap(cur_ind_colors,:);
    
    set(cScatter, 'CData', cur_colors);
    title(sprintf('Time Step: %d', time_pt*3));
    drawnow;
    pause(0.18);

end


%% 3.1 REAL DATA - Create 3D movie/animation of scatterplots

% event = 6;

trial_batch_analysis = 0; % TODO change this manually (= 1 if analyzing multiple inspection events, =0 if individual)
alignments = {'outcome','inspection','response','trialstart'};
%window_nums = {1,2}; % 1 = pre-alignment, 2 = post-alignment (e.g., pre-outcome = anticipation, post-outcome = feedback)

subject_IDs = {'EMU024','EMU030'};
for subject_num = 1:2
    subject_ID = subject_IDs{subject_num};
   
    final_out_dir = sprintf('/media/Data/Human_Intracranial_MAD/analysis/Cross_Freq_Coupling/%s/%s/synchrony_data/generalized_phase',subject_ID, reference);
    out_dir = sprintf('/media/Data/Human_Intracranial_MAD/analysis/Cross_Freq_Coupling/%s/%s/synchrony_data/generalized_phase',subject_ID, reference);
    if ~exist(out_dir, 'dir'), mkdir(out_dir); end
    
    load(sprintf('%s/%s/%s_MAD_SES%d_Setup.mat',data_base_dir,subject_ID,subject_ID,sesnum));
    
    % load brainstorm file and outer contacts
    cTable = readtable(sprintf('/media/Data/Human_Intracranial_MAD/analysis/Location_Consistency/3D Coordinates/%s.tsv',subject_ID), "FileType","delimitedtext",'Delimiter', '\t');
    
    channel_names_bs = cTable.Channel; % extract "Channels" variable in table
    coord3d_all = cTable.SCS;
    
    channel_names_list = [];
    for i = 1:numel(channel_names_bs)
        cur_name = channel_names_bs{i};
        channel_names_list = [channel_names_list; convertCharsToStrings(cur_name)];
    end
    
    all_contact_coord = [];
    for i = 1:numel(elec_name)
        idx = find(channel_names_list == convertCharsToStrings(elec_name{i})); %outer contacts values are strings
        if isempty(idx), continue; end
        coord3d = coord3d_all{idx};
    
        coord3d(1) = ''; coord3d(numel(coord3d)) = '';
        commas = find(coord3d == ',');
        x_coord = coord3d(1:commas(1)-1); x_coord = str2double(convertCharsToStrings(x_coord));
        y_coord = coord3d(commas(1)+1:commas(2)-1); y_coord = str2double(convertCharsToStrings(y_coord));
        z_coord = coord3d(commas(2)+1:numel(coord3d)); z_coord = str2double(convertCharsToStrings(z_coord));
        all_contact_coord = [all_contact_coord; [x_coord y_coord z_coord]];
    
    end
    
    for align_iter = 2 % CHANGE THIS MANUALLY
        alignment = alignments{align_iter};
    
        switch alignment
            case 'outcome'
                windows = {'anticipation','feedback'};
                window_num = 2; %TODO: manually change number
                window = windows{window_num}; % feedback
                time_window_seconds = 1; % set length for aligned data window
                time_window_length = time_window_seconds*Fs; % set sampling rate-dependent spectral parameters
                if window_num == 1
                    win = [time_window_length 0];
                else
                    win = [0 time_window_length];
                end
            case 'inspection'
                % only post-alignment period matters
                %window = 'post-inspection';
                
                %[time_btwn_inspections, last_inspection_to_response] = filter_inspections(subject_ID, sesnum);
    
                time_window_seconds = 0.5; % TODO: get_time_between_inspections.m
                time_window_length = time_window_seconds*Fs; % set sampling rate-dependent spectral parameters
                win = [0 time_window_length];
            case 'response'
                windows = {'pre-response','post-response'};
                window = windows{1}; % TODO: manually change number
                % TODO: time_window_length differs depending on pre- or post-alignment
                time_window_seconds = 0.5;
                time_window_length = time_window_seconds*Fs; % set sampling rate-dependent spectral parameters
                if window_num == 1
                    win = [time_window_length 0];
                else
                    win = [0 time_window_length];
                end
            case 'trialstart' 
                window = 'pre-trial'; % to analyze inter-trial intervals!
                time_window_seconds = 0.5;
                time_window_length = time_window_seconds*Fs;
                win = [time_window_length 0];
        end
    
        addpath('/media/Data/Human_Intracranial_MAD/_toolbox/spectral-analysis');
        [align_times,trial_numbers] = get_align_times(filters, trial_times, trial_words, alignment);
    
        remove_ind = isnan(align_times);
        align_times(isnan(align_times)) = [];
        trial_numbers(remove_ind) = [];
        
        align_times = round(align_times*Fs); % convert alignment times from seconds to samples
    
        % some random events to plot
        for event = [1 2 3 4 5 10 31] % TEMP 5/22/25 added
    
            xgp_full_res = [];
            xgp_full_extra = [];
            xgp_all_contacts = []; % each column = 1 (outer) contact, each row = 1 timept
    
            for contact_iter = 1:numel(elec_name)
        
                if ~ismember(elec_name(contact_iter), channel_names_bs)
                    continue;
                end
                
                load(sprintf('%s/extracted_phase/xgp_%s_ses%02d_ch%03d.mat',final_out_dir,elec_name{contact_iter},sesnum,contact_iter));
                % variable loaded: xgp
                
                if trial_batch_analysis == 1
                    trial_num = 3; % TODO change this manually
                    align_times_trial = align_times(trial_numbers == trial_num); % get all align_times (inspections) associated with this trial number
                    last_idx_trial = find(trial_numbers == trial_num, 1, "last");
                    time_window_samples = align_times(last_idx_trial + 1) - align_times_trial(1);
                    win = [0 time_window_samples];
                    xgp_window = xgp((round(align_times_trial(1)-win(1)):round(align_times_trial(1)+win(2)-1)));
                    xgp_to_plot = xgp_window(3:3:numel(xgp_window));
                
                else
                    %event = 2; % TODO change this manually             % TEMP
                    %for event = 1:numel(align_times) 
                    
                    % if (strcmp(alignment, 'response') && strcmp(window, 'pre-response'))
                    %     time_window_samples = last_inspection_to_response(event);
                    %     win = [time_window_samples 0]; 
                    % elseif strcmp(alignment, 'inspection') && trial_numbers(event + 1) ~= trial_numbers(event)
                    %     % if this event is the last inspection in a trial, use time between last inspection and response
                    %     time_window_samples = last_inspection_to_response(trial_numbers(event));
                    %     win = [0 time_window_samples];
                    % elseif strcmp(alignment, 'inspection')
                    %     time_window_samples = time_btwn_inspections(trial_numbers(event));
                    %     win = [0 time_window_samples]; % post-inspection
                    % elseif strcmp(alignment, 'trialstart')
                    %     time_window_seconds = 0.5;
                    %     time_window_length = time_window_seconds*Fs;
                    %     win = [time_window_length 0];
                    % end
            
                    xgp_window = xgp((round(align_times(event)-win(1)):round(align_times(event)+win(2)-1)));
                    xgp_window_extra = xgp((round(align_times(event)-win(1)):round(align_times(event)+2*win(2))));
                    xgp_to_plot = xgp_window(3:3:numel(xgp_window)); % start:time_step:num samples in time window = sampling rate * num seconds)
            
                end
        
                xgp_full_res = [xgp_full_res xgp_window];
                xgp_full_extra = [xgp_full_extra xgp_window_extra];
                xgp_all_contacts = [xgp_all_contacts xgp_to_plot]; % HORIZONTAL
        
            end % end loop through all contacts
    
            clear angle; % clear the variable called "angle"
            angles_full_extra = angle(xgp_full_extra);
            angles_full_res = angles_full_extra(1:height(xgp_full_res),:);
            angles_all_contacts = angle(xgp_all_contacts);
        
            % save the data to plot
            fname_to_plot = sprintf('/media/Data/Human_Intracranial_MAD/analysis/Cross_Freq_Coupling/%s/Ground/synchrony_data/generalized_phase/data_to_plot/%s_event%03d.mat',subject_ID,alignment, event);
            save(fname_to_plot, "angles_all_contacts")
        
        end % TEMP 5/22/25 added
    
    end
end % end loop through subjects
% clear angle; % clear the variable called "angle"
% angles_all_contacts = angle(xgp_all_contacts);
% samples_to_zero = zeros(numel(xgp_all_contacts),1); % TODO


% load data to plot from this directory:
    % TEMP!!
alignment = 'inspection'; % trialstart % alignment to load
event = 5; % event to load
data_to_plot_dir = '/media/Data/Human_Intracranial_MAD/analysis/Cross_Freq_Coupling/EMU025/Ground/synchrony_data/generalized_phase/data_to_plot';
load(sprintf('%s/%s_event%03d.mat',data_to_plot_dir,alignment,event))
%if strcmp(subject_ID, 'EMU025'), angles_all_contacts = angles_all_contacts(:,1:206); end

% 1. loop through time samples to create matrix of c_indices (repeat for additional alignment points)
    % normalize each ROW of angles_all_contacts
addpath('/media/Data/Human_Intracranial_MAD/analysis/PAC_code/CF_Coupling/generalized phase');
%cmap = colorcet('C5'); % colorcet('C2','N',500); % 500 different colors
cmap = blue_cyclic_colormap(256);
% 'C2' is cyclic rainbow, 'C5' is monochrome
num_colors = size(cmap, 1);

colors_ind_all = zeros(height(angles_all_contacts),width(angles_all_contacts)); % each column = 1 outer contact, each row = 1 timept

for i = 1:height(angles_all_contacts) % iterate through time points
    
    cur_timept_angles = angles_all_contacts(i,:);

    ph_min = min(cur_timept_angles); ph_max = max(cur_timept_angles);
    c_norm_angles = (cur_timept_angles - ph_min) / (ph_max - ph_min); % convert from [-pi pi] to [0 1] % TODO: check that this works!!
    c_indices = round(c_norm_angles * (num_colors - 1)) + 1; % convert from [0 1] to [0 256], since 256 different colors in cmap
    colors_ind_all(i,:) = c_indices; % use these indices to access the cmap variable
    
    %c_angles = cmap(c_indices,:); % these are the COLORS that will be used as argument in scatter() TODO: rename this to c_colors
    %colors_all = [colors_all; c_angles]; 

end

% 2. loop through time samples again to generate colored scatter plots BASED ON matrix colors_all
% plot row-by-row of angles_all_contacts/colors_all - each row = 1 frame
% create movie here

% initialize first frame
figure('Position', [300, 50, 750, 750]);
cur_ind_colors = colors_ind_all(1,:); 
cur_colors = cmap(cur_ind_colors,:);
cScatter = scatter3(all_contact_coord(:,1), all_contact_coord(:,2), all_contact_coord(:,3), 100, cur_colors, 'filled'); % 200 is marker size

colormap(cmap) % Apply cyclic colormap
colorbar
clim([-pi pi]) % Set colorbar range to match phase angles % changed from caxis to clim

if strcmp(subject_ID,'EMU024')
    xlim([-30, 45])
    ylim([-75, 10])
    zlim([5, 75])
    view(15, 40);
elseif strcmp(subject_ID,'EMU036')
    xlim([-45, 45])
    ylim([-5, 80])
    zlim([20, 90])
    view(15, 40);  % set viewing angle: Azimuth (rotation around z-axis) = 15°, Elevation (angle above xy plane) = 40°
elseif strcmp(subject_ID,'EMU030')
    xlim([-30, 40])
    ylim([-5, 85])
    zlim([15, 90])
    view(20, 40);
end

xlabel('X') 
ylabel('Y') 
zlabel('Z')

% rest of the frames
for time_pt = 2:height(angles_all_contacts)
    cur_ind_colors = colors_ind_all(time_pt,:); %cur_angles = angles_all_contacts(time_pt,:);
    cur_colors = cmap(cur_ind_colors,:);
    
    set(cScatter, 'CData', cur_colors);
    title(sprintf('Time Step: %d', time_pt*3));
    drawnow;
    pause(0.1);

end

%% 3.2 Visualizing wavefront - 2D, binary colors

wavefront = zeros(size(angles_all_contacts));
for i = 1:height(angles_all_contacts)
    wavefront(i,:) = (abs(angles_all_contacts(i,:) - 3.14)) < 0.15;
    %wavefront(i,:) = (angles_all_contacts(i,:) == min(angles_all_contacts(i,:)));
end

cur_colors = zeros(15,3);
cur_wavefront = wavefront(1,:);
for i = 1:numel(cur_wavefront)
    if cur_wavefront(i) == 1
        cur_colors(i,:) = [1 0 0];
    else
        cur_colors(i,:) = [0 0 1];
    end
end

figure;
%cScatter = scatter(outer_contact_coord(:,1), 30*ones(size(outer_contact_coord(:,3))), 200, cur_colors, 'filled');
cScatter = scatter(outer_contact_coord(:,1), outer_contact_coord(:,3), 200, cur_colors, 'filled');
xlim([-45, 45]) 
ylim([20, 90])
xlabel('X') 
ylabel('Z')

for time_pt = 2:height(angles_all_contacts)
    cur_wavefront = wavefront(time_pt,:);
    for i = 1:numel(cur_wavefront)
        if cur_wavefront(i) == 1
            cur_colors(i,:) = [1 0 0];
        else
            cur_colors(i,:) = [0 0 1];
        end
    end
    
    set(cScatter, 'CData', cur_colors);
    title(sprintf('Time Step: %d', time_pt*3));
    drawnow;
    pause(0.18);

end


%% TEMP: reconcile channel_name contacts with cTable names

reference = 'neighbor_average';
subject_ID = 'EMU024';
sesnum = 1;

addpath("/media/Data/Human_Intracranial_MAD/_toolbox/spectral-analysis")

datapath = '/media/Data/Human_Intracranial_MAD/';
data_base_dir = sprintf('%s1_formatted',datapath);
out_dir = sprintf('/media/Data/Human_Intracranial_MAD/analysis/Cross_Freq_Coupling/%s',subject_ID);
if ~exist(out_dir, 'dir'), mkdir(out_dir); end

% load setup files
load(sprintf('%s/%s/%s_MAD_SES%d_Setup_%s.mat',data_base_dir,subject_ID,subject_ID,sesnum,reference));

% load all contact 3d coordinates
cTable = readtable(sprintf('/media/Data/Human_Intracranial_MAD/analysis/Location_Consistency/3D Coordinates/%s.tsv',subject_ID), "FileType","delimitedtext",'Delimiter', '\t');
channel_names = cTable.Channel;

channel_nameNOTbrainstorm = setdiff(channel_name, channel_names);
brainstormNOTchannel_name = setdiff(channel_names, channel_name);

if ~isempty(channel_nameNOTbrainstorm)
    save(sprintf('%s/channel_nameNOTINbrainstorm',out_dir),"channel_nameNOTbrainstorm");
end
if ~isempty(brainstormNOTchannel_name)
    save(sprintf('%s/brainstormNOTINchannel_name',out_dir),"brainstormNOTchannel_name");
end

% END TEMP SECTION
