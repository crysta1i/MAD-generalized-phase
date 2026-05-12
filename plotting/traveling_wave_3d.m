% Animation of 3d traveling wave

subject_ID = 'EMU025';
reference = 'Ground';
sesnum = 1;

datapath = '/media/Data/Human_Intracranial_MAD/';
data_base_dir = sprintf('%s1_formatted',datapath);
final_out_dir = sprintf('/media/Data/Human_Intracranial_MAD/analysis/Cross_Freq_Coupling/%s/%s/synchrony_data/generalized_phase',subject_ID, reference);
out_dir = sprintf('/media/Data/Human_Intracranial_MAD/analysis/Cross_Freq_Coupling/%s/%s/synchrony_data/generalized_phase',subject_ID, reference);

load(sprintf('%s/%s/%s_MAD_SES%d_Setup.mat',data_base_dir,subject_ID,subject_ID,sesnum));

cTable = readtable(sprintf('/media/Data/Human_Intracranial_MAD/analysis/Location_Consistency/3D Coordinates/%s.tsv',subject_ID), "FileType","delimitedtext",'Delimiter', '\t');
channel_names_bs = cTable.Channel;
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
if strcmp(subject_ID, 'EMU025'), all_contact_coord = all_contact_coord(1:206,:); end

alignment = 'inspection'; % trialstart % alignment to load
event = 5; % event to load [3 5 9]
data_to_plot_dir = '/media/Data/Human_Intracranial_MAD/analysis/Cross_Freq_Coupling/EMU025/Ground/synchrony_data/generalized_phase/data_to_plot';
load(sprintf('%s/%s_event%03d.mat',data_to_plot_dir,alignment,event))
if strcmp(subject_ID, 'EMU025'), angles_all_contacts = angles_all_contacts(:,1:206); end


addpath('/media/Data/Human_Intracranial_MAD/analysis/PAC_code/CF_Coupling/generalized phase');
addpath('/media/Data/Human_Intracranial_MAD/analysis/PAC_code/CF_Coupling/generalized phase/circstat-matlab');
cmap = blue_cyclic_colormap(256);
num_colors = size(cmap, 1);

colors_ind_all = zeros(height(angles_all_contacts),width(angles_all_contacts)); % each column = 1 outer contact, each row = 1 timept

for i = 1:height(angles_all_contacts) % iterate through time points
    
    cur_timept_angles = angles_all_contacts(i,:);

    ph_min = min(cur_timept_angles); ph_max = max(cur_timept_angles);
    c_norm_angles = (cur_timept_angles - ph_min) / (ph_max - ph_min); % convert from [-pi pi] to [0 1] % TODO: check that this works!!
    c_indices = round(c_norm_angles * (num_colors - 1)) + 1; % convert from [0 1] to [0 256], since 256 different colors in cmap
    colors_ind_all(i,:) = c_indices; 

end

% initialize first frame
figure('Position', [300, 50, 750, 750]);
cur_ind_colors = colors_ind_all(1,:); 
cur_colors = cmap(cur_ind_colors,:);
cScatter = scatter3(all_contact_coord(:,1), all_contact_coord(:,2), all_contact_coord(:,3), 100, cur_colors, 'filled'); % 200 is marker size

colormap(cmap)
colorbar
clim([-pi pi])

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
    cur_ind_colors = colors_ind_all(time_pt,:);
    cur_colors = cmap(cur_ind_colors,:);
    
    set(cScatter, 'CData', cur_colors);
    title(sprintf('Time Step: %d', time_pt*3));
    drawnow;
    pause(0.1);

end