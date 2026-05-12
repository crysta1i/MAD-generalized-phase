% gp for MAD study (in PAC_code/CF_coupling/generalized phase)
% save extracted generalized phase, code reused from gp_demo.m

% 5/27/25: 5-40Hz data extracted for EMU036, 030, 024
%          5-20Hz data extracted for EMU036

brain_regions = ["amygdala","hippocampus","accumbens","orbitofrontal", "insula","anterior cingulate","middle cingulate","posterior cingulate","caudate","dorsolateral frontal", "mesial frontal", "dorsolateral parietal","dorsolateral temporal", "mesial temporal", "temporal pole", "dorsal occipital", "mesial occipital"];
brain_region_names = ["amygdala","hippocampus","accumbens","orbitofrontal", "insula","cingulate_ant","cingulate_mid","cingulate_post","caudate","dors_frontal", "mes_frontal", "dors_parietal","dors_temp", "mes_temp", "temp_pole", "dors_occ", "mes_occ"];

reference = 'Ground'; 
subject_IDs = {'EMU001','EMU024', 'EMU025','EMU030', 'EMU036', 'EMU037','EMU038','EMU039','EMU041','EMU051'};
for subject_num = 3 %[2 3 4 5 8] % can't run on subj 1 because have not saved ind_regions_all, etc
    %TODO: EMU041 doesn't work, 51 and 37 don't have brainstorm data (can't determine outer contact)
    % EMU038 (subjnum 7) is problematic because of different naming ->
    % channel_ind creation couldn't find appropriate indices
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
    
    %paths
    datapath = '/media/Data/Human_Intracranial_MAD/';
    data_base_dir = sprintf('%s1_formatted',datapath);

    % Create generalized phase directory
    out_dir = sprintf('/media/Data/Human_Intracranial_MAD/analysis/Cross_Freq_Coupling/%s/%s/synchrony_data',subject_ID, reference);
    final_out_dir = sprintf('/media/Data/Human_Intracranial_MAD/analysis/Cross_Freq_Coupling/%s/%s/synchrony_data/generalized_phase',subject_ID, reference);
    if ~exist(out_dir, 'dir'), mkdir(out_dir); end

    for sesnum = 1 %:num_sessions %TEMP
        % ONLY WORKS if reference = Ground
        if strcmp(subject_ID, 'EMU036') && sesnum == 2
                continue;
        else
            %load(sprintf('%s/%s/%s_MAD_SES%d_Setup_%s.mat',data_base_dir,subject_ID,subject_ID,sesnum,reference));
            load(sprintf('%s/%s/%s_MAD_SES%d_Setup.mat',data_base_dir,subject_ID,subject_ID,sesnum));
            %if strcmp(subject_ID, 'EMU025'), elec_name = elec_name(1:206); end
            % channel_name_na = channel_name;
            % channel_location_na = channel_location;
            % channel_ind_na = channel_ind;
        end

        %% Compute gp on only brain areas / channels of interest
        
        % 1.  2/20/25: Outermost contacts on each electrode

        % Load ind_regions_all 
        fname = sprintf('/media/Data/Human_Intracranial_MAD/analysis/Cross_Freq_Coupling/%s/%s/synchrony_data/index_regions_relevantChannels.mat',subject_ID, 'neighbor_average');
        load(fname); load(sprintf('/media/Data/Human_Intracranial_MAD/analysis/Cross_Freq_Coupling/%s/%s/synchrony_data/regions_relevantChannels.mat',subject_ID, 'neighbor_average'));
        
        %[outer_contacts, outer_contacts_ind] = get_outer_contact(subject_ID,reference);
        load(sprintf('%s/extracted_phase/outer_contacts/outer_contact_names.mat',final_out_dir)); % outer_contacts
        load(sprintf('%s/extracted_phase/outer_contacts/outer_contacts_ind.mat',final_out_dir)); % outer_contacts_ind
       
        for i = 1:numel(outer_contacts_ind) % 1:numel(elec_name)
            
            
            ind_cur = outer_contacts_ind(i); % TEMP
            
            % find where ind_regions_all == outer_contacts_ind(i), and get corresponding regions_all element
            %ind_idx = find(ind_regions_all == outer_contacts_ind(i)); %TEMP
            %region_name = regions_all(ind_idx); %TEMP
 
            % Only works for GROUND reference
            % TEMP commented out:
            %load(sprintf('%s/%s/separate_channel_files/%s_MAD_SES%d_ch%03d.mat',data_base_dir,subject_ID,subject_ID,sesnum,ind_cur)); % for saving phase for outer contacts only
            load(sprintf('%s/%s/separate_channel_files/%s_MAD_SES%d_ch%03d.mat',data_base_dir,subject_ID,subject_ID,sesnum,i)); % for saving phase for ALL contacts
            
            %load(sprintf('%s/%s/separate_channel_files/%s/%s_MAD_SES%d_ch%03d.mat',data_base_dir,subject_ID,reference,subject_ID,sesnum,ind_cur)); % replaced chnum_na with ind_cur
                
            %% FILTERING

            addpath('/media/Data/Human_Intracranial_MAD/analysis/PAC_code/CF_Coupling/generalized phase/generalized-phase/analysis');

            % parameters
            filter_order = 4; Fs = 2048; lp = 5;
            dt = 1 / Fs; T = length(data) / Fs; time = dt:dt:T;
            
            % preprocessing - lowpass + notch
            [b,a] = butter( 4, [5 200] ./ (Fs/2) ); xw = filtfilt( b, a, data ); 
            [b,a] = butter( 4, [58 62] ./ (Fs/2), 'stop' ); xw = filtfilt( b, a, xw ); 
            [b,a] = butter( 4, [115 125] ./ (Fs/2), 'stop' ); xw = filtfilt( b, a, xw );
            
            % create GP
            [b,a] = butter( filter_order, [5 40] ./ (Fs/2) ); xf = filtfilt( b, a, data ); 
            xgp = generalized_phase_vector( xf, Fs, lp );

            % SAVE
            % if ~exist(sprintf('%s/extracted_phase/outer_contacts',final_out_dir),"dir")
            %     mkdir(sprintf('%s/extracted_phase/outer_contacts',final_out_dir));
            % end
            fname = sprintf('%s/extracted_phase/outer_contacts/xgp_ses%02d_ch%03d.mat',final_out_dir,sesnum,ind_cur);
            %fname = sprintf('%s/extracted_phase/xgp_%s_ses%02d_ch%03d.mat',final_out_dir,elec_name{i},sesnum,i);
            
            save(fname,"xgp")

            % NEW - DID NOT RUN YET:
            %fname = sprintf('%s/extracted_phase/outer_contacts/xw_%s_ses%02d_ch%03d.mat',final_out_dir,region_name,sesnum,ind_cur);
            %save(fname,"xw")

            %fname = sprintf('%s/extracted_phase/outer_contacts/xf_%s_ses%02d_ch%03d.mat',final_out_dir,region_name,sesnum,ind_cur);
            %save(fname,"xf")
            

        end % loop through outer_contacts_ind

    end % loop through sessions

end % loop through subjects






%% plot - wideband LFP and GP

fg1 = figure; hold on; ax1 = gca; set( fg1, 'position', [ 88  1593  1250  420 ] )
plot( time, xw, 'linewidth', 4, 'color', 'k' ); h4 = cline( time, xf, [], angle(xgp) );
set( h4, 'linestyle', '-', 'linewidth', 5  ); xlim( [0.06 .65] ); axis off
l1 = line( [.1 .2], [-125 -125] ); set( l1, 'linewidth', 4, 'color', 'k' )
l2 = line( [.1 .1], [-125 -75] ); set( l2, 'linewidth', 4, 'color', 'k' )



%%

fg1 = figure; hold on; ax1 = gca; set( fg1, 'position', [ 88  1593  1250  420 ] )
for i = 1:2

    load(sprintf('%s/%s/separate_channel_files/%s_MAD_SES%d_ch%03d.mat', ...
        data_base_dir,subject_ID,subject_ID,sesnum,cur_elec_contact_ind(i))); % i = contact number on current electrode
    
    filter_order = 4; Fs = 2048; lp = 5;
    dt = 1 / Fs; T = length(data) / Fs; time = dt:dt:T;
    
    [b,a] = butter( filter_order, [5 40] ./ (Fs/2) ); xf = filtfilt( b, a, data ); 
    
    [b,a] = butter( 4, [5 200] ./ (Fs/2) ); xw1 = filtfilt( b, a, data ); 
    [b,a] = butter( 4, [58 62] ./ (Fs/2), 'stop' ); xw1 = filtfilt( b, a, xw1 ); 
    [b,a] = butter( 4, [115 125] ./ (Fs/2), 'stop' ); xw1 = filtfilt( b, a, xw1 );
    
    % main figure
    plot( time, xw1, 'linewidth', 2, 'color', 'k' ); 
    h4 = cline( time, xf, [], angle(xgp_cur_elec(:,1)) );
    set( h4, 'linestyle', '-', 'linewidth', 3  ); 

end

xw1 = xw1(round(align_times(event)-win(1)):round(align_times(event)+win(2)-1));

%xlim( [0.06 .65] ); axis off
%l1 = line( [.1 .2], [-125 -125] ); set( l1, 'linewidth', 4, 'color', 'k' )
%l2 = line( [.1 .1], [-125 -75] ); set( l2, 'linewidth', 4, 'color', 'k' )

% inset
% map = colorcet( 'C2' ); colormap( circshift( map, [ 28, 0 ] ) )
% ax2 = axes; set( ax2, 'position', [0.2116    0.6976    0.0884    0.2000] ); axis image
% [x1,y1] = pol2cart( angle( exp(1i.*linspace(-pi,pi,100)) ), ones( 1, 100 ) );
% h3 = cline( x1, y1, linspace(-pi,pi,100) ); axis off; set( h3, 'linewidth', 6 )

% text labels
% t1 = text( 0, 0, 'GP' );
% set( t1, 'fontname', 'arial', 'fontsize', 28, 'fontweight', 'bold', 'horizontalalignment', 'center' )
% set( gcf, 'currentaxes', ax1 )
% t2 = text( 0.1260, -146.8832, '100 ms' );
% set( t2, 'fontname', 'arial', 'fontsize', 24, 'fontweight', 'bold' )
% t2 = text( 0.0852, -130.2651, '50 \muV' );
% set( t2, 'fontname', 'arial', 'fontsize', 24, 'fontweight', 'bold', 'rotation', 90 )

for subject_num = [2 3 4 5 6 7 8 9 10] %TODO: EMU041 doesn't work (rerun later)
    subject_ID = subject_IDs{subject_num};
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

    datapath = '/media/Data/Human_Intracranial_MAD/';
    data_base_dir = sprintf('%s1_formatted',datapath);

    % Create generalized phase directory
    out_dir = sprintf('/media/Data/Human_Intracranial_MAD/analysis/Cross_Freq_Coupling/%s/%s/synchrony_data',subject_ID, reference);
    final_out_dir = sprintf('/media/Data/Human_Intracranial_MAD/analysis/Cross_Freq_Coupling/%s/%s/synchrony_data/generalized_phase',subject_ID, reference);
    if ~exist(out_dir, 'dir'), mkdir(out_dir); end

    [outer_contacts, outer_contacts_ind] = get_outer_contact(subject_ID,reference);

    if ~exist(sprintf('%s/extracted_phase/outer_contacts',final_out_dir),"dir")
        mkdir(sprintf('%s/extracted_phase/outer_contacts',final_out_dir));
    end
    fname_contacts = sprintf('%s/extracted_phase/outer_contacts/outer_contact_names.mat',final_out_dir);
    fname_contact_ind = sprintf('%s/extracted_phase/outer_contacts/outer_contacts_ind.mat',final_out_dir);
    save(fname_contacts,"outer_contacts")
    save(fname_contact_ind,"outer_contacts_ind")
    
end

t = readtable("/media/Data/Human_Intracranial_MAD/analysis/Location_Consisteny/3D Coordinates/EMU036.tsv", "FileType","delimitedtext",'Delimiter', '\t');

%% Fix file naming (5/20/25)

% Set folder
folder = '/media/Data/Human_Intracranial_MAD/analysis/Cross_Freq_Coupling/EMU036/Ground/synchrony_data/generalized_phase/extracted_phase';

% Get file list
files = dir(fullfile(folder, '*.*'));
files = files(~[files.isdir]);

% Loop through each file
for k = 1:length(files)
    oldName = files(k).name;
    [~, name, ext] = fileparts(oldName);
    
    newDigits = '';
    matched = false;

    for i = 1:numel(elec_name)
        if i <= 153, continue; end
        
        if contains(name, sprintf('_%s_',elec_name{i}))
            newDigits = i;
            matched = true;
            break;
        end
    end
    
    % Skip file if no matching code was found
    if ~matched
        fprintf('No matching code for: %s\n', oldName);
        continue;
    end
    
    % Replace the last 3 digits with the new value
    newName = strrep(name, '153', sprintf('%03d',newDigits));  % replace final 3 digits

    % Full paths
    oldFullPath = fullfile(folder, oldName);
    newFullPath = fullfile(folder, [newName ext]);

    % Rename file
    movefile(oldFullPath, newFullPath);
end

disp('Renaming complete!');
