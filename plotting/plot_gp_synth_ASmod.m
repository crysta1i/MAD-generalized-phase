%% 0. Synthetic Data

%%% AS modification %%%
% It seems like you only use this to define a number of points, so I'm
% converting it to seconds with an assumed Fs of 100 Hz. (Just for
% convenience below.)
Fs = 100;
time = (1:300)/Fs;
clear angle;
%%%%%%%%%%%%%%%%%%%%%%%

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


%%% AS modifications %%%
% Data setup: Wavefront moving
    % horizontal: sin { (x + (0 - min_x)) * (7/ (max_x - min_x) ) + (0.1 * t) }
spatial_freq = 1; % spatial frequency of the wave in periods per unit distance
propagation_speed = 1; % propagation speed of wave in periods per second
wavenumber = 2*pi*spatial_freq; % converted to radians per unit distance
prop_speed_rad = wavenumber*propagation_speed; % propagation speed in radians per second

synth_horiz = zeros(numel(time),height(outer_contact_coord));
for i = 1:numel(time)
    for j = 1:height(outer_contact_coord)
        x_coor = outer_contact_coord(j,1);
        norm_x = (x_coor + 45) / 90; % max_x - min_x = 90, % 0 - min_x = 45
        t = time(i);
        synth_horiz(i, j) = sin(wavenumber*norm_x - prop_speed_rad*t);
        % synth_horiz(i, j) = sin(2*pi*norm_x + ((3*2*pi)/300)*t);
        %synth_horiz(i, j) = sin( (x_coor + 45)*( 2*pi /(90)) + 0.3*i );
    end
end
phase_horiz = asin(synth_horiz); 
data = phase_horiz;
%%%%%%%%%%%%%%%%%%%%%%%%

    % vertical: sin { (y - min_y) * (7/ (max_y - min_y) ) + (0.1 * t) } % t = i
synth_vert = zeros(numel(time),height(outer_contact_coord));
for i = 1:numel(time)
    for j = 1:height(outer_contact_coord)
        z_coor = outer_contact_coord(j,3); % 3rd column
        synth_vert(i, j) = sin((z_coor - 20)*(7/(70)) + 0.3*i); % TODO
    end
end
phase_vert = asin(synth_vert);
data = phase_vert;


%% ALL DATA: create visualization:

addpath('/media/Data/Human_Intracranial_MAD/analysis/PAC_code/CF_Coupling/generalized phase');

blue_cyclic_colormap = [[zeros(256,1) zeros(256,1) linspace(0,1,256)']; [[linspace(0,0.5,128)'; linspace(0.5,0,128)'] [linspace(0,0.5,128)'; linspace(0.5,0,128)'] linspace(1,0,256)']];
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
    cur_ind_colors = colors_ind_all(time_pt,:); %cur_angles = angles_all_contacts(time_pt,:);
    cur_colors = cmap(cur_ind_colors,:);
    
    set(cScatter, 'CData', cur_colors);
    title(sprintf('Time Step: %d', time_pt));
    drawnow;
    pause(0.18);

end