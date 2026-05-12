% project contact coordinates (2d) onto line

% Prereqs: 
% 1. outer_contact_coord (produced in plot_gp.m section 1)
% 2. angles_all_contacts (produced in plot_gp.m section 3) -- OUTER CONTACTS (2D VERSION) of angles_all_contacts
% 3. align_times -- inspection events and inter-trial events extracted (produced in plot_gp.m)

%% CHOOSE A LINE
% Find circle and radius around the data points

% range of X/Y/Z
rangeX = max(outer_contact_coord(:,1)) - min(outer_contact_coord(:,1));
rangeY = max(outer_contact_coord(:,2)) - min(outer_contact_coord(:,2));
rangeZ = max(outer_contact_coord(:,3)) - min(outer_contact_coord(:,3));

% get radius
radiusXY = 0.6 * max([rangeX rangeY]); % if only using x and y coordinates
radiusXZ = 0.6 * max([rangeX rangeZ]); % if only using x and z coordinates


centerX = 0.5 * (max(outer_contact_coord(:,1)) + min(outer_contact_coord(:,1)));
centerY = 0.5 * (max(outer_contact_coord(:,2)) + min(outer_contact_coord(:,2))); 
centerZ = 0.5 * (max(outer_contact_coord(:,3)) + min(outer_contact_coord(:,3)));

%% 0.1 Quantifying traveling wave 
% OUTER CONTACTS
% prereqs for running this section: outer_contacts, outer_contact_coord, outer_contacts_ind, 
% setup file, align_times, win, event, channel_names_bs, final_out_dir,sesnum 

xgp_all_contacts = [];
for contact_iter = 1:numel(outer_contacts)
    if ~ismember(outer_contacts(contact_iter), channel_names_bs), continue; end

    ind_idx = outer_contacts_ind(contact_iter);
    
    load(sprintf('%s/extracted_phase/xgp_%s_ses%02d_ch%03d.mat',final_out_dir,elec_name{ind_idx},sesnum,ind_idx));
    % variable loaded: xgp

    xgp_window = xgp((round(align_times(event)-win(1)):round(align_times(event)+win(2)-1)));
    xgp_to_plot = xgp_window(3:3:numel(xgp_window));

    xgp_all_contacts = [xgp_all_contacts xgp_to_plot];

end % end loop through outer contacts

angles_all_contacts = angle(xgp_all_contacts);


angles = [0, 1/8, 1/4, 3/8, 1/2, 5/8, 3/4, 7/8]; % scale pi by these amounts
angle_names = {'0','1/8','1/4', '3/8', '1/2', '5/8', '3/4', '7/8'};

figure; % 8 subplots (1 created per iteration of below loop) in 1 figure
tiledlayout(4,2);

for i = 1:numel(angles)
    theta = angles(i);
    theta_name = angle_names{i};

    all_dists_on_proj = zeros(height(angles_all_contacts),1); % distance on projection, at each time pt
    num_time_points = height(angles_all_contacts); % angles_all_contacts is already a subset of all timepoints (e.g., every 3 points)
    for j = 1:num_time_points
    
        % find contact that is closest to 0 phase
        [~, wavefront_contact] = min(abs(angles_all_contacts(j,:)));

    
        lineStart_theta = (theta) * pi;
        lineEnd_theta = lineStart_theta + pi;
        
        % convert (theta, R) form polar to cartesian
        lineStartX = centerX + radiusXZ * cos(lineStart_theta);
        lineStartY = centerZ + radiusXZ * sin(lineStart_theta);
        lineEndX = centerX + radiusXZ * cos(lineEnd_theta);
        lineEndY = centerZ + radiusXZ * sin(lineEnd_theta);
        
        lineStart = [lineStartX, lineStartY];
        lineEnd = [lineEndX, lineEndY];
        
        % PROJECT POINTS ONTO THE LINE
        
        % if only using X and Z
        pointX = outer_contact_coord(wavefront_contact,1);
        pointY = outer_contact_coord(wavefront_contact,3);
        point = [pointX, pointY];
    
        vec_line = lineEnd - lineStart;
        point_to_line = point - lineStart;
        projection = lineStart + ( (point_to_line * vec_line') / (vec_line * vec_line') ) * vec_line;
    
        % get dist from projection to start of line
        x_dist = projection(1) - lineStart(1);
        y_dist = projection(2) - lineStart(2);
        dist = sqrt(x_dist.^2 + y_dist.^2);
            
        all_dists_on_proj(j) = dist;        
    
    end % end loop through time points

    nexttile;
    plot(all_dists_on_proj,'.') % make 8 subplots in one figure 
    title(sprintf("%s pi", theta_name))

end % end loop through angles

%% 0.2 Quantifying traveling wave 
% ALL CONTACTS -- INDIV ELECTRODES

% track where the wavefront is located on each electrode across time
% Prereqs: setup, alignment, win, cTable (of coordinates) loaded

elec_letters = [];
for i = 1:numel(channel_names_bs)
    cur_name = channel_names_bs{i};
    if contains(cur_name, "'")
        if ~ismember(cur_name(1:2),elec_letters)
            elec_letters = [elec_letters; convertCharsToStrings(cur_name(1:2))];
        end
    else
        if ~ismember(cur_name(1),elec_letters)
            elec_letters = [elec_letters; convertCharsToStrings(cur_name(1))];
        end
    end
end
   
for event = [1 2 3 4 5 10 29 38 101 123 399 420]
    figure;
    tlplot = tiledlayout(3,3);
    for i = 1:9 %numel(elec_letters)
        cur_letter = elec_letters(i);
        xgp_cur_elec = []; % each column = 1 contact on current electrode, each row = 1 timept
        xgp_elec_to_plot = [];
        j = 1;
        if ~contains(cur_letter, "'")
            while ~contains(elec_name{j},cur_letter) || contains(elec_name{j},"'"), j = j+1; end
        else
            while ~contains(elec_name{j},cur_letter), j = j+1; end
        end
        cur_name = elec_name{j};
        cur_el    % load setup fileec_contact_names = []; % stores name of contacts on current electrode
        cur_elec_contact_ind = []; % stores indices of contacts on curr elec
        while strcmp(cur_name(1:2),cur_letter) || strcmp(cur_name(1),cur_letter)
            cur_elec_contact_names = [cur_elec_contact_names; convertCharsToStrings(cur_name)];
            cur_elec_contact_ind = [cur_elec_contact_ind; j];
            load(sprintf('%s/extracted_phase/xgp_%s_ses%02d_ch%03d.mat',final_out_dir,cur_name,sesnum,j));
            
            xgp_cur_contact = xgp((round(align_times(event)-win(1)):round(align_times(event)+win(2)-1)));
            xgp_cur_contact_to_plot = xgp_cur_contact(3:3:numel(xgp_cur_contact));
            
            xgp_cur_elec = [xgp_cur_elec xgp_cur_contact];
            xgp_elec_to_plot = [xgp_elec_to_plot xgp_cur_contact_to_plot];
            
            j = j+1; 
            if j > numel(elec_name), break; end
            cur_name = elec_name{j};
        end
    
        angles_cur_elec = angle(xgp_cur_elec);
        angles_elec_to_plot = angles_cur_elec(3:3:height(angles_cur_elec),:);
    

        % select smaller time windows to analyze
        subwindow = [860 960];
        angles_cur_elec = angles_cur_elec(subwindow(1):subwindow(2),:);
        peak_times = zeros(1,width(angles_cur_elec));
        for e = 1:width(angles_cur_elec)
            [peak1, peak_time] = min(abs(angles_cur_elec(:,e)-pi));
            % if peak_time == 1 || peak_time == numel(subwindow(1):subwindow(2))
            %     peak_time_est = peak_time;
            % elseif angles_cur_elec(peak_time,e) * angles_cur_elec(peak_time + 1,e) <= 0
            %     peak_time_est = interp1([angles_cur_elec(peak_time,e), angles_cur_elec(peak_time + 1,e)], [peak_time, peak_time + 1], 0, 'linear');
            % else
            %     peak_time_est = interp1([angles_cur_elec(peak_time - 1,e), angles_cur_elec(peak_time,e)], [peak_time - 1, peak_time], 0, 'linear');
            % end
            % peak_times(e) = peak_time_est;
            peak_times(e) = peak_time;
        end
            %peak_times = peak_times - min(peak_times);
        figure;
        plot(peak_times,'.','MarkerSize',25);
        xlabel("contact num")
        ylabel("time samples to trough")
        title("Electrode X (0.43-0.48s)")


        % get coordinates of the contacts in current electrode
            % TEMP
        coord3d_all = cTable.SCS;
        coord_cur_elec = [];
        delete_columns = [];
        for k = 1:numel(cur_elec_contact_names)
            idx = find(channel_names_list == convertCharsToStrings(cur_elec_contact_names(k))); %outer contacts values are strings
            if isempty(idx)
                delete_columns = [delete_columns k];
                continue; 
            end
            coord3d = coord3d_all{idx};
    
            coord3d(1) = ''; coord3d(numel(coord3d)) = '';
            commas = find(coord3d == ',');
            x_coord = coord3d(1:commas(1)-1); x_coord = str2double(convertCharsToStrings(x_coord));
            y_coord = coord3d(commas(1)+1:commas(2)-1); y_coord = str2double(convertCharsToStrings(y_coord));
            z_coord = coord3d(commas(2)+1:numel(coord3d)); z_coord = str2double(convertCharsToStrings(z_coord));
            coord_cur_elec = [coord_cur_elec; [x_coord y_coord z_coord]];
    
        end
        angles_elec_to_plot(:,delete_columns) = [];
    
        % distance of wavefront on current electrode, at each time pt
        all_dists = zeros(height(angles_elec_to_plot),1); 
        num_time_points = height(angles_elec_to_plot);
        for jj = 1:num_time_points
        
            % find contact that is closest to 0 phase
            [~, wavefront_contact] = min(abs(angles_elec_to_plot(jj,:)));
    
            % find distance between "last" coordinate on current electrode and wavefront
            x_dist = coord_cur_elec(wavefront_contact,1) - coord_cur_elec(1,1);
            y_dist = coord_cur_elec(wavefront_contact,2) - coord_cur_elec(1,2);
            z_dist = coord_cur_elec(wavefront_contact,3) - coord_cur_elec(1,3);
            dist = sqrt(x_dist.^2 + y_dist.^2 + z_dist.^2);
            all_dists(jj) = dist;
    
        end
    
        nexttile;
        plot(all_dists,'.')
        title(sprintf('%s',cur_letter))

        % correlation between phase and distance
        corrs_sig_all_times = zeros(height(angles_elec_to_plot),2);
        for kk = 1:height(angles_elec_to_plot)
            pl = angles_elec_to_plot(kk,:); % j is time point index
            [ cc, pv ] = circ_corrcl(pl, 1:numel(cur_elec_contact_names)); % all in a row so exact coordinates don't matter
            corrs_sig_all_times(kk,1) = cc;
            corrs_sig_all_times(kk,2) = pv < 0.05;
        end

        % longest streak of significant timepoints
        maxsum = 0; locmax = 1; sum = 0;
        for l = 1:numel(corrs_sig_all_times(:,2))
            if corrs_sig_all_times(l,2) == 1
                sum = sum + 1; 
            else
                if sum > maxsum, maxsum = sum; locmax = l; end
                sum = 0;
            end
        end
        longest_streak = maxsum; % index/timepoint right after this streak is locmax
    end

    title(tlplot, 'Wavefront Location vs Time')
    tlplot.TileSpacing = 'compact';
    tlplot.Padding = 'compact';

    fname = sprintf('%s/single_elec_plots/%s_ses%02d_event%03d.jpg',final_out_dir,alignment,sesnum,event);
    if ~exist(sprintf('%s/single_elec_plots',final_out_dir),'dir')
        mkdir(sprintf('%s/single_elec_plots',final_out_dir));
    end
    print(gcf,'-djpeg', fname)
    close;
end

%% Plot neighboring contacts' signals to visualize phase offsets

% plot(angles_elec_to_plot(1:100,1))
% hold on
% plot(angles_elec_to_plot(1:100,2))
% plot(angles_elec_to_plot(1:100,3)) % 3 at a time
% hold off
alignment = 'trialstart'; %inspection,trialstart
win = [1024 0]; %[0 1024] [1024 0]

[align_times,trial_numbers] = get_align_times(filters, trial_times, trial_words, alignment);
remove_ind = isnan(align_times);
align_times(isnan(align_times)) = [];
trial_numbers(remove_ind) = [];
align_times = round(align_times*Fs);

% event = 5;

plot_out_dir = sprintf('/media/Data/Human_Intracranial_MAD/analysis/Cross_Freq_Coupling/%s/Ground/synchrony_data/generalized_phase/waveform_plots',subject_ID);
if ~exist(plot_out_dir,'dir'), mkdir(plot_out_dir); end

for event = [1 3 4 5 10 12]
    for i = 1:numel(elec_letters)
        
        cur_letter = elec_letters(i);
        xgp_cur_elec = []; % each column = 1 contact on current electrode, each row = 1 timept
        xgp_elec_to_plot = [];
        j = 1;
        if ~contains(cur_letter, "'")
            while ~contains(elec_name{j},cur_letter) || contains(elec_name{j},"'"), j = j+1; end
        else
            while ~contains(elec_name{j},cur_letter), j = j+1; end
        end
        cur_name = elec_name{j};
        cur_elec_contact_names = []; % stores name of contacts on current electrode
        cur_elec_contact_ind = []; % stores indices of contacts on curr elec
        while strcmp(cur_name(1:2),cur_letter) || strcmp(cur_name(1),cur_letter)
            cur_elec_contact_names = [cur_elec_contact_names; convertCharsToStrings(cur_name)];
            cur_elec_contact_ind = [cur_elec_contact_ind; j];
            load(sprintf('%s/extracted_phase/xgp_%s_ses%02d_ch%03d.mat',final_out_dir,cur_name,sesnum,j));
            
            xgp_cur_contact = xgp((round(align_times(event)-win(1)):round(align_times(event)+win(2)-1)));
            xgp_cur_contact_to_plot = xgp_cur_contact(3:3:numel(xgp_cur_contact));
            
            xgp_cur_elec = [xgp_cur_elec xgp_cur_contact];
            xgp_elec_to_plot = [xgp_elec_to_plot xgp_cur_contact_to_plot];
            
            j = j+1; 
            if j > numel(elec_name), break; end
            cur_name = elec_name{j};
        end
    
        angles_cur_elec = angle(xgp_cur_elec);
        angles_elec_to_plot = angles_cur_elec(3:3:height(angles_cur_elec),:);
    
        fg1 = figure; hold on; ax1 = gca; set( fg1, 'position', [ 88  1593  1250  420 ] )
        for contact_iter = 1:numel(cur_elec_contact_ind)
        
            load(sprintf('%s/%s/separate_channel_files/%s_MAD_SES%d_ch%03d.mat', ...
                data_base_dir,subject_ID,subject_ID,sesnum,cur_elec_contact_ind(contact_iter))); 
            
            filter_order = 4; Fs = 2048; lp = 5;
            dt = 1 / Fs; T = length(data) / Fs; time = dt:dt:T;
            
            [b,a] = butter( filter_order, [5 40] ./ (Fs/2) ); xf = filtfilt( b, a, data ); 
            
            [b,a] = butter( 4, [5 200] ./ (Fs/2) ); xw1 = filtfilt( b, a, data ); 
            [b,a] = butter( 4, [58 62] ./ (Fs/2), 'stop' ); xw1 = filtfilt( b, a, xw1 ); 
            [b,a] = butter( 4, [115 125] ./ (Fs/2), 'stop' ); xw1 = filtfilt( b, a, xw1 );
            
            xf = xf(round(align_times(event)-win(1)):round(align_times(event)+win(2)-1));
            time = time(round(align_times(event)-win(1)):round(align_times(event)+win(2)-1));
            time = time - time(1);
            if numel(cur_elec_contact_names) < 11
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
                plot(time, xf, 'linewidth', 3, 'color', colors(contact_iter,:));
            else
                colors = [
                    0.80 0.90 1.00;
                    0.76 0.86 0.98;
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
                    0.06 0.26 0.63;
                    0.00 0.22 0.60];
                plot(time, xf, 'linewidth', 3, 'color', colors(contact_iter+(19-numel(cur_elec_contact_names)),:));
            end 
        end
        title(sprintf('Electrode %s',cur_letter))
        fname = sprintf('%s/elec%s_%s_event%03d.jpg', plot_out_dir, cur_letter, alignment, event);
        print(gcf,'-djpeg', fname)
        close;
    end
end


fg1 = figure; hold on; ax1 = gca; set( fg1, 'position', [ 88  1593  1250  420 ] )
for contact_iter = 1:15

    load(sprintf('%s/%s/separate_channel_files/%s_MAD_SES%d_ch%03d.mat', ...
        data_base_dir,subject_ID,subject_ID,sesnum,cur_elec_contact_ind(contact_iter))); 
    
    filter_order = 4; Fs = 2048; lp = 5;
    dt = 1 / Fs; T = length(data) / Fs; time = dt:dt:T;
    
    [b,a] = butter( filter_order, [5 40] ./ (Fs/2) ); xf = filtfilt( b, a, data ); 
    
    [b,a] = butter( 4, [5 200] ./ (Fs/2) ); xw1 = filtfilt( b, a, data ); 
    [b,a] = butter( 4, [58 62] ./ (Fs/2), 'stop' ); xw1 = filtfilt( b, a, xw1 ); 
    [b,a] = butter( 4, [115 125] ./ (Fs/2), 'stop' ); xw1 = filtfilt( b, a, xw1 );
    
    xf = xf(round(align_times(event)-win(1)):round(align_times(event)+win(2)-1));
    time = time(round(align_times(event)-win(1)):round(align_times(event)+win(2)-1));
    time = time - time(1);
    % main figure
    %colors = ['k','b','r','g','y','m','c'];
    % colors = [
    %     0.80 0.90 1.00;
    %     0.71 0.81 0.96;
    %     0.63 0.73 0.92;
    %     0.54 0.64 0.87;
    %     0.45 0.56 0.83;
    %     0.36 0.47 0.78;
    %     0.27 0.39 0.74;
    %     0.18 0.30 0.69;
    %     0.09 0.22 0.65;
    %     0.00 0.13 0.60];
    colors = [
        0.80 0.90 1.00;
        0.76 0.86 0.98;
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
        0.06 0.26 0.63;
        0.00 0.22 0.60];
    %plot( time, xw1, 'linewidth', 2, 'color', colors(i) ); 
    % h4 = cline( time, xf, [], angle(xgp_cur_elec(:,1)) );
    % set( h4, 'linestyle', '-', 'linewidth', 3  ); 
    plot(time, xf, 'linewidth', 3, 'color', colors(contact_iter+3,:)); % 
end


%% 1. Visualize results for one event at a time (dot plots across time)

% Statistics to store
% 1. visualize results across ANGLES (histogram: how many timepoints had each angle be significant)
% 2. visualize results across TIME PTS (scatterplot: for each time point, plot the angle that had highest corrcoef)
    % do something special to the points that are significant

angles = [0, 1/8, 1/4, 3/8, 1/2, 5/8, 3/4, 7/8]; % scale pi by these amounts
    
xgp_all_contacts = [];

for contact_iter = 1:numel(outer_contacts)

    %ind_cur = outer_contacts_ind(contact_iter); % TEMP (for outer contacts)

    % find where ind_regions_all == outer_contacts_ind(i)
        % TEMP commented out (for outer contacts)
    %ind_idx = find(ind_regions_all == outer_contacts_ind(contact_iter));
    %region_name = regions_all(ind_idx);
    %load(sprintf('%s/extracted_phase/outer_contacts/xgp_%s_ses%02d_ch%03d.mat',final_out_dir,region_name,sesnum,ind_cur));
    
    if ~ismember(outer_contacts(contact_iter), channel_names_bs)
        continue;
    end

    ind_idx = outer_contacts_ind(contact_iter);
    
    load(sprintf('%s/extracted_phase/xgp_%s_ses%02d_ch%03d.mat',final_out_dir,elec_name{ind_idx},sesnum,ind_idx));
    % variable loaded: xgp



    % % define WIN here, using times btwn inspections (and redo all dot plots??)
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
    % end



    xgp_window = xgp((round(align_times(event)-win(1)):round(align_times(event)+win(2)-1)));
    xgp_to_plot = xgp_window(3:3:numel(xgp_window));

    xgp_all_contacts = [xgp_all_contacts xgp_to_plot];

end % end loop through outer contacts

angles_all_contacts = angle(xgp_all_contacts);

% 2. project contacts, get statistics 
significant_angles = zeros(1, numel(angles)); % same size as angles
max_corr_angles = []; % angle with max correlation at each time point
max_corr_angles_pv = []; % tracks whether above angles are significant

    % n_points = 1000;  % Number of points to plot the circle smoothly
    % theta = linspace(0, 2*pi, n_points);   % Angle from 0 to 2π
    % x = centerX + radiusXZ * cos(theta);   % x = x0 + r*cos(θ)
    % y = centerZ + radiusXZ * sin(theta);
    % 
    % scatter(outer_contact_coord(:,1), outer_contact_coord(:,3), 200);
    % hold on;
    % plot(x, y, 'b-', 'LineWidth', 2);
    % hold off;

num_time_points = height(angles_all_contacts); % angles_all_contacts is already a subset of all timepoints (e.g., every 3 points)
for j = 1:num_time_points

    cur_timept_max_corr = -1;
    cur_timept_max_corr_angle = -1;
    cur_timept_max_corr_pv = -1;

    % choose a theta angle
    for i = 1:numel(angles)
        theta = angles(i);
        lineStart_theta = (theta) * pi;
        lineEnd_theta = lineStart_theta + pi;
        
        % convert (theta, R) form polar to cartesian
        lineStartX = centerX + radiusXZ * cos(lineStart_theta);
        lineStartY = centerZ + radiusXZ * sin(lineStart_theta);
        lineEndX = centerX + radiusXZ * cos(lineEnd_theta);
        lineEndY = centerZ + radiusXZ * sin(lineEnd_theta);
        
        lineStart = [lineStartX, lineStartY];
        lineEnd = [lineEndX, lineEndY];
        
        %% PROJECT POINTS ONTO THE LINE
        
        dist_on_projection = zeros(height(outer_contact_coord),1);
        for k = 1:height(outer_contact_coord)
            % if only using X and Z
            pointX = outer_contact_coord(k,1);
            pointY = outer_contact_coord(k,3);
            point = [pointX, pointY];
        
            vec_line = lineEnd - lineStart;
            point_to_line = point - lineStart;
            projection = lineStart + ( (point_to_line * vec_line') / (vec_line * vec_line') ) * vec_line;
        
            % get dist from projection to start of line
            x_dist = projection(1) - lineStart(1);
            y_dist = projection(2) - lineStart(2);
            dist = sqrt(x_dist.^2 + y_dist.^2);
            dist_on_projection(k) = dist;
        end
        
        %% CORRELATION BETWEEN DISTANCE AND PHASE
        pl = angles_all_contacts(j,:); % j is time point index
        dist_on_projection( isnan(pl) ) = []; pl( isnan(pl) ) = [];
        [ cc, pv ] = circ_corrcl(pl, dist_on_projection);
        if pv < 0.05
            significant_angles(i) = significant_angles(i) + 1; % increment the count
        end
        if cc > cur_timept_max_corr
            cur_timept_max_corr = cc;
            cur_timept_max_corr_angle = theta;
            cur_timept_max_corr_pv = pv;
        end
    end % end loop through angles
    max_corr_angles = [max_corr_angles; cur_timept_max_corr_angle];
    max_corr_angles_pv = [max_corr_angles_pv; cur_timept_max_corr_pv];
end % end loop through time points
max_corr_angles_sig = max_corr_angles_pv < 0.05;
corr_pv = [max_corr_angles max_corr_angles_sig];

% 3. visualize results
x = 1:numel(max_corr_angles);
plot(x(max_corr_angles_sig), max_corr_angles(max_corr_angles_sig),'r.')
hold on
plot(x(~max_corr_angles_sig), max_corr_angles(~max_corr_angles_sig),'b.')
xlabel('Time point');
ylabel('Angle of max correlation');
legend('Red = significant', 'Blue = ns');
hold off;

%histogram(max_corr_angles(max_corr_angles_sig)) 

%% 2. Visualize results for multiple events (histogram)

% Go through all sub-alignments (inspection types) in loop

% copied from plot_gp script: %

subject_ID = 'EMU036'; sesnum = 1; reference = 'Ground'; Fs = 2048;
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

% load ind_regions_all and regions_all
% NOTE: here, reference = neighbor_average because these variables were saved in neighbor_avg folder
fname = sprintf('/media/Data/Human_Intracranial_MAD/analysis/Cross_Freq_Coupling/%s/%s/synchrony_data/index_regions_relevantChannels.mat',subject_ID, 'neighbor_average');
load(fname);
load(sprintf('/media/Data/Human_Intracranial_MAD/analysis/Cross_Freq_Coupling/%s/%s/synchrony_data/regions_relevantChannels.mat',subject_ID, 'neighbor_average'));   

% load setup files
load(sprintf('%s/%s/%s_MAD_SES%d_Setup.mat',data_base_dir,subject_ID,subject_ID,sesnum));

% load xgp for each contact of interest, extract data for one event on each xgp data vector, concatenate data for all contacts

% 0. loop through outer contacts
xgp_all_contacts = []; % each column = 1 (outer) contact, each row = 1 timept

align_iter = 2; % CHANGE THIS MANUALLY
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

if strcmp(alignment, 'inspection')
    inspect_types = {'win_domain','loss_domain','amount','probability','novel_attribute_inspection','repeated_inspection',...
                'fourth_unique_attribute','first_unique_attribute'};
    % 'first_unique_attribute','second_unique_attribute','third_unique_attribute',
   
   
    for align2_iter = 1:numel(inspect_types)
        alignment2 = inspect_types{align2_iter};

        [align_times,trial_numbers] = get_align_times(filters, trial_times, trial_words, alignment2);
        remove_ind = isnan(align_times);
        align_times(isnan(align_times)) = [];
        trial_numbers(remove_ind) = [];

        align_times = round(align_times*Fs); % convert alignment times from seconds to samples
        
        angles = [0, 1/8, 1/4, 3/8, 1/2, 5/8, 3/4, 7/8]; % scale pi by these amounts
        
        significant_angles = zeros(1, numel(angles)); 
        max_corr_angles = []; % angle with maximum correlation at each time point
        max_corr_angles_pv = []; % tracks whether above angles are significant
        event_nums = 1:50; %randi(numel(align_times),50,1); % 10 random events
        %TEMP: NON-random event selection 
        
        for event_iter = 1:numel(event_nums)
            
            event = event_nums(event_iter);
            xgp_all_contacts = []; % one for each event
        
            % 1. get appropriate angles_all_contacts for the given event#
            for contact_iter = 1:numel(outer_contacts)

                if ~ismember(outer_contacts(contact_iter), channel_names_bs)
                    continue;
                end
        
                ind_idx = outer_contacts_ind(contact_iter);
                
                load(sprintf('%s/extracted_phase/xgp_%s_ses%02d_ch%03d.mat',final_out_dir,elec_name{ind_idx},sesnum,ind_idx));
                % variable loaded: xgp
        
                [time_btwn_inspections, last_inspection_to_response] = filter_inspections(subject_ID, sesnum);
        
                % TODO: DEFINE WIN HERE AGAIN!! using times btwn inspections (and redo all histograms)
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
                % end

                time_window_samples = 0.5 * Fs;
                win = [0 time_window_samples];

                xgp_window = xgp((round(align_times(event)-win(1)):round(align_times(event)+win(2)-1)));
                xgp_to_plot = xgp_window(3:3:numel(xgp_window));
        
                xgp_all_contacts = [xgp_all_contacts xgp_to_plot];
        
            end % end loop through outer contacts
        
            angles_all_contacts = angle(xgp_all_contacts);
        
            % 2. project contacts, get statistics 
            
                % n_points = 1000;  % Number of points to plot the circle smoothly
                % theta = linspace(0, 2*pi, n_points);   % Angle from 0 to 2π
                % x = centerX + radiusXZ * cos(theta);   % x = x0 + r*cos(θ)
                % y = centerZ + radiusXZ * sin(theta);
                % 
                % scatter(outer_contact_coord(:,1), outer_contact_coord(:,3), 200);
                % hold on;
                % plot(x, y, 'b-', 'LineWidth', 2);
                % hold off;
            
            num_time_points = height(angles_all_contacts); % angles_all_contacts is already a subset of all timepoints (e.g., every 3 points)
            for j = 1:num_time_points
            
                cur_timept_max_corr = -1;
                cur_timept_max_corr_angle = -1;
                cur_timept_max_corr_pv = -1;
            
                % choose a theta angle
                for i = 1:numel(angles)
                    theta = angles(i);
                    lineStart_theta = (theta) * pi;
                    lineEnd_theta = lineStart_theta + pi;
                    
                    % convert (theta, R) form polar to cartesian
                    lineStartX = centerX + radiusXZ * cos(lineStart_theta);
                    lineStartY = centerZ + radiusXZ * sin(lineStart_theta);
                    lineEndX = centerX + radiusXZ * cos(lineEnd_theta);
                    lineEndY = centerZ + radiusXZ * sin(lineEnd_theta);
                    
                    lineStart = [lineStartX, lineStartY];
                    lineEnd = [lineEndX, lineEndY];
                    
                    %% PROJECT POINTS ONTO THE LINE
                    
                    dist_on_projection = zeros(height(outer_contact_coord),1);
                    for k = 1:height(outer_contact_coord)
                        % if only using X and Z
                        pointX = outer_contact_coord(k,1);
                        pointY = outer_contact_coord(k,3);
                        point = [pointX, pointY];
                    
                        vec_line = lineEnd - lineStart;
                        point_to_line = point - lineStart;
                        projection = lineStart + ( (point_to_line * vec_line') / (vec_line * vec_line') ) * vec_line;
                    
                        % get dist from projection to start of line
                        x_dist = projection(1) - lineStart(1);
                        y_dist = projection(2) - lineStart(2);
                        dist = sqrt(x_dist.^2 + y_dist.^2);
                        dist_on_projection(k) = dist;
                    end
                    
                    %% CORRELATION BETWEEN DISTANCE AND PHASE
                    pl = angles_all_contacts(j,:); % j is time point index
                    dist_on_projection( isnan(pl) ) = []; pl( isnan(pl) ) = [];
                    [ cc, pv ] = circ_corrcl(pl, dist_on_projection);
                    if pv < 0.05
                        significant_angles(i) = significant_angles(i) + 1; % increment the count
                    end
                    if cc > cur_timept_max_corr
                        cur_timept_max_corr = cc;
                        cur_timept_max_corr_angle = theta;
                        cur_timept_max_corr_pv = pv;
                    end
                end % end loop through angles
                max_corr_angles = [max_corr_angles; cur_timept_max_corr_angle];
                max_corr_angles_pv = [max_corr_angles_pv; cur_timept_max_corr_pv];
            end % end loop through time points

        
        end
        
        %corr_pv = [max_corr_angles max_corr_angles_sig];
        max_corr_angles_sig = max_corr_angles_pv < 0.05;
        
        figure;
        histogram(max_corr_angles(max_corr_angles_sig))
        title(sprintf('%s',alignment2))
        fname = sprintf('/media/Data/Human_Intracranial_MAD/analysis/Cross_Freq_Coupling/EMU036/Ground/synchrony_data/generalized_phase/dist_phase_corr_histograms/%s_halfsec.jpg', ...
            alignment2);
        print(gcf,'-djpeg', fname)
        close;

    end

end



% copied from plot_gp script^ %


%% One sub-alignment at a time

angles = [0, 1/8, 1/4, 3/8, 1/2, 5/8, 3/4, 7/8]; % scale pi by these amounts

% f_inspect = figure;
% xlabel('Time point');
% ylabel('Angle of max correlation');
% legend('Red = significant', 'Blue = ns');

significant_angles = zeros(1, numel(angles)); % preallocate; same size as angles
max_corr_angles = []; % angle with maximum correlation at each time point
max_corr_angles_pv = []; % tracks whether above angles are significant
event_nums = 1:randi(numel(align_times),50,1); % 10 random events
% put everything below into a loop through events %

for event_iter = 1:numel(event_nums)
    
    event = event_nums(event_iter);
    xgp_all_contacts = []; % one for each event

    % 1. get appropriate angles_all_contacts for the given event#
    for contact_iter = 1:numel(outer_contacts)

        %ind_cur = outer_contacts_ind(contact_iter); % TEMP (for outer contacts)
    
        % find where ind_regions_all == outer_contacts_ind(i)
            % TEMP commented out (for outer contacts)
        %ind_idx = find(ind_regions_all == outer_contacts_ind(contact_iter));
        %region_name = regions_all(ind_idx);
        %load(sprintf('%s/extracted_phase/outer_contacts/xgp_%s_ses%02d_ch%03d.mat',final_out_dir,region_name,sesnum,ind_cur));
        
        if ~ismember(outer_contacts(contact_iter), channel_names_bs)
            continue;
        end

        ind_idx = outer_contacts_ind(contact_iter);
        
        load(sprintf('%s/extracted_phase/xgp_%s_ses%02d_ch%03d.mat',final_out_dir,elec_name{ind_idx},sesnum,ind_idx));
        % variable loaded: xgp

        [time_btwn_inspections, last_inspection_to_response] = filter_inspections(subject_ID, sesnum);

        % TODO: DEFINE WIN HERE AGAIN!! using times btwn inspections (and redo all histograms)
        if (strcmp(alignment, 'response') && strcmp(window, 'pre-response'))
            time_window_samples = last_inspection_to_response(event);
            win = [time_window_samples 0]; 
        elseif strcmp(alignment, 'inspection') && trial_numbers(event + 1) ~= trial_numbers(event)
            % if this event is the last inspection in a trial, use time between last inspection and response
            time_window_samples = last_inspection_to_response(trial_numbers(event));
            win = [0 time_window_samples];
        elseif strcmp(alignment, 'inspection')
            time_window_samples = time_btwn_inspections(trial_numbers(event));
            win = [0 time_window_samples]; % post-inspection
        end


        xgp_window = xgp((round(align_times(event)-win(1)):round(align_times(event)+win(2)-1)));
        xgp_to_plot = xgp_window(3:3:numel(xgp_window));

        xgp_all_contacts = [xgp_all_contacts xgp_to_plot];

    end % end loop through outer contacts

    angles_all_contacts = angle(xgp_all_contacts);

    % 2. project contacts, get statistics 
    
        % n_points = 1000;  % Number of points to plot the circle smoothly
        % theta = linspace(0, 2*pi, n_points);   % Angle from 0 to 2π
        % x = centerX + radiusXZ * cos(theta);   % x = x0 + r*cos(θ)
        % y = centerZ + radiusXZ * sin(theta);
        % 
        % scatter(outer_contact_coord(:,1), outer_contact_coord(:,3), 200);
        % hold on;
        % plot(x, y, 'b-', 'LineWidth', 2);
        % hold off;
    
    num_time_points = height(angles_all_contacts); % angles_all_contacts is already a subset of all timepoints (e.g., every 3 points)
    for j = 1:num_time_points
    
        cur_timept_max_corr = -1;
        cur_timept_max_corr_angle = -1;
        cur_timept_max_corr_pv = -1;
    
        % choose a theta angle
        for i = 1:numel(angles)
            theta = angles(i);
            lineStart_theta = (theta) * pi;
            lineEnd_theta = lineStart_theta + pi;
            
            % convert (theta, R) form polar to cartesian
            lineStartX = centerX + radiusXZ * cos(lineStart_theta);
            lineStartY = centerZ + radiusXZ * sin(lineStart_theta);
            lineEndX = centerX + radiusXZ * cos(lineEnd_theta);
            lineEndY = centerZ + radiusXZ * sin(lineEnd_theta);
            
            lineStart = [lineStartX, lineStartY];
            lineEnd = [lineEndX, lineEndY];
            
            %% PROJECT POINTS ONTO THE LINE
            
            dist_on_projection = zeros(height(outer_contact_coord),1);
            for k = 1:height(outer_contact_coord)
                % if only using X and Z
                pointX = outer_contact_coord(k,1);
                pointY = outer_contact_coord(k,3);
                point = [pointX, pointY];
            
                vec_line = lineEnd - lineStart;
                point_to_line = point - lineStart;
                projection = lineStart + ( (point_to_line * vec_line') / (vec_line * vec_line') ) * vec_line;
            
                % get dist from projection to start of line
                x_dist = projection(1) - lineStart(1);
                y_dist = projection(2) - lineStart(2);
                dist = sqrt(x_dist.^2 + y_dist.^2);
                dist_on_projection(k) = dist;
            end
            
            %% CORRELATION BETWEEN DISTANCE AND PHASE
            pl = angles_all_contacts(j,:); % j is time point index
            dist_on_projection( isnan(pl) ) = []; pl( isnan(pl) ) = [];
            [ cc, pv ] = circ_corrcl(pl, dist_on_projection);
            if pv < 0.05
                significant_angles(i) = significant_angles(i) + 1; % increment the count
            end
            if cc > cur_timept_max_corr
                cur_timept_max_corr = cc;
                cur_timept_max_corr_angle = theta;
                cur_timept_max_corr_pv = pv;
            end
        end % end loop through angles
        max_corr_angles = [max_corr_angles; cur_timept_max_corr_angle];
        max_corr_angles_pv = [max_corr_angles_pv; cur_timept_max_corr_pv];
    end % end loop through time points
    
    
    % 3. visualize results
    
    % x = 1:numel(max_corr_angles);
    % plot(x(max_corr_angles_sig), max_corr_angles(max_corr_angles_sig),'r.')
    % hold on
    % plot(x(~max_corr_angles_sig), max_corr_angles(~max_corr_angles_sig),'b.')
    % hold off;

end

%corr_pv = [max_corr_angles max_corr_angles_sig];
max_corr_angles_sig = max_corr_angles_pv < 0.05;

histogram(max_corr_angles(max_corr_angles_sig)) 

%% 3. Determine whether wavefront is traveling
% Plot location of wavefront (outer contacts)

% RUN PLOT_GP 3.0 (until 2nd breakpt) AND SECTION 1 OF THIS SCRIPT
% BEFORE RUNNING THIS SECTION (need max_corr_angles_sig variable)

% in each row of angles_all_contacts, IF phase-distance was significant at that time point, 
% find the contact(s) that were closest to 0 phase

contacts_to_plot = zeros(height(angles_all_contacts),2); % column1 = x-coord, column2 = z-coord
for i = 1:height(angles_all_contacts)
    [~, idx_contact] = min(abs(angles_all_contacts(i,:)));
    contacts_to_plot(i,:) = [outer_contact_coord(idx_contact,1), outer_contact_coord(idx_contact,3)];
end

times = 1:height(angles_all_contacts);
x_to_plot = contacts_to_plot(:,1);
y_to_plot = contacts_to_plot(:,2);

% IF significant at that time point, plot in red. O/W, plot in blue.

% time vs x and z coordinate of wavefront
plot3(times(max_corr_angles_sig), x_to_plot(max_corr_angles_sig), y_to_plot(max_corr_angles_sig), 'r.')
hold on
plot3(times(~max_corr_angles_sig), x_to_plot(~max_corr_angles_sig), y_to_plot(~max_corr_angles_sig), 'b.')
ylim([-45, 45])
zlim([20, 90])
xlabel('Time point');
ylabel('x coord');
zlabel('z coord');
%legend('Red = significant', 'Blue = ns');
hold off;

% time vs x-coordinate of wavefront
figure;
plot(times(max_corr_angles_sig), x_to_plot(max_corr_angles_sig), 'r.')
hold on
plot(times(~max_corr_angles_sig), x_to_plot(~max_corr_angles_sig), 'b.')
ylim([-45, 45])
xlabel('Time point');
ylabel('x coord');
%legend('Red = significant', 'Blue = ns');
hold off;

% time vs z-coordinate of wavefront
figure;
plot(times(max_corr_angles_sig), y_to_plot(max_corr_angles_sig), 'r.')
hold on
plot(times(~max_corr_angles_sig), y_to_plot(~max_corr_angles_sig), 'b.')
ylim([20, 90])
xlabel('Time point');
ylabel('z coord');
%legend('Red = significant', 'Blue = ns');
hold off;
