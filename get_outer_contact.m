function [outer_contacts, outer_contacts_ind] = get_outer_contact(subject_ID, reference)
% Find outermost contact on each electrode
% outer_contacts: names of contacts
% outer_contacts_ind: corresponding indices in setup files

% if not enough args
if nargin < 2
    reference = 'Ground';
end

sesnum = 1;

datapath = '/media/Data/Human_Intracranial_MAD/';
data_base_dir = sprintf('%s1_formatted',datapath);

%% TODO: Loop through subjects

subject_IDs = {'EMU001','EMU024', 'EMU025','EMU030', 'EMU036', 'EMU037','EMU038','EMU039','EMU041','EMU047','EMU051'};

for subject_num = [1 2 3 4 5 8 9] % 7 is EMU038
    subject_ID = subject_IDs{subject_num};

    % find all unique electrode letters in the subject (loop through elec_name)
    if strcmp(reference, 'Ground')
        load(sprintf('%s/%s/%s_MAD_SES%d_Setup.mat',data_base_dir,subject_ID,subject_ID,1),"elec_name","elec_ind","elec_area"); % set sesnum = 1
    else
        load(sprintf('%s/%s/%s_MAD_SES%d_Setup_%s.mat',data_base_dir,subject_ID,subject_ID,1,reference),"channel_name","channel_ind","channel_area");
    end
    
    cTable = readtable(sprintf('/media/Data/Human_Intracranial_MAD/analysis/Location_Consistency/3D Coordinates/%s.tsv',subject_ID), "FileType","delimitedtext",'Delimiter', '\t');
    channel_names_bs = cTable.Channel; % bs = brainstorm
    
    elec_letters = [];
    %elec_nums = [];
    
    for i = 1:numel(elec_name)
    
        cur_name = elec_name{i};
        
        if contains(cur_name,"'")
            cur_elec_char = cur_name(1:2);
        else
            cur_elec_char = cur_name(1);
        end
    
        elec_letters = [elec_letters; convertCharsToStrings(cur_elec_char)];
    
    end
    elec_letters = unique(elec_letters);
    % now we have list of elec_letters
    

    outer_contacts = [];
    %outer_contacts_su = []; outer_contacts_ind_su = [];
    %outer_contacts_bs = []; % do not have or update an outer_contacts_bs_ind var
    for i = 1:numel(elec_letters)
    
        % extract all items in elec_name
        names_cur_elec_su = []; %names_cur_elec_ind_su = [];
        nums_cur_elec_su = []; % elec names but omit the letters (just the nums)
        %ind_cur_elec_su = []; % corresponding indices
        for j = 1:numel(elec_name)
            cur_name = elec_name{j};

            if strcmp(elec_area{j},"Out")
                continue;
            end

            if contains(elec_letters(i),"'")

                if contains(cur_name,elec_letters(i))
                    names_cur_elec_su = [names_cur_elec_su; convertCharsToStrings(elec_name{j})];
                    %ind_cur_elec_su = [ind_cur_elec_su; elec_ind(j)];
                    
                    if contains(cur_name,"'")
                        if strcmp(subject_ID, 'EMU038')
                            cur_num = cur_name(7:numel(cur_name));
                        else
                            cur_num = cur_name(3:numel(cur_name));
                        end
                    else
                        cur_num = cur_name(2:numel(cur_name));
                    end
        
                    nums_cur_elec_su = [nums_cur_elec_su; str2double(convertCharsToStrings(cur_num))];
                end

            else

                if contains(cur_name,elec_letters(i)) && ~contains(cur_name,"'")
                    names_cur_elec_su = [names_cur_elec_su; convertCharsToStrings(elec_name{j})];
                    %names_cur_elec_ind_su = [names_cur_elec_ind_su; elec_ind(j)];
                    
                    % if contains(cur_name,"'")
                    %     if strcmp(subject_ID, 'EMU038')
                    %         cur_num = cur_name(7:numel(cur_name));
                    %     else
                    %         cur_num = cur_name(3:numel(cur_name));
                    %     end
                    % else
                    cur_num = cur_name(2:numel(cur_name));
                    % end
        
                    nums_cur_elec_su = [nums_cur_elec_su; str2double(convertCharsToStrings(cur_num))];
                end

            end
        end
        [cur_outer_su, idx_cur_outer_su] = max(nums_cur_elec_su); % max for THIS elec letter, setup
        cur_outer_name_su = names_cur_elec_su(idx_cur_outer_su);
        %outer_contacts_su = [outer_contacts_su; cur_outer_su];
        %outer_contacts_ind_su = [outer_contacts_ind_su; idx_cur_outer_su];
    
        % extract all items in channel_names_bs
        names_cur_elec_bs = [];
        nums_cur_elec_bs = [];
        for j = 1:numel(channel_names_bs)
            cur_name = channel_names_bs{j};

            if contains(elec_letters(i),"'")
                
                if contains(cur_name,elec_letters(i))
    
                    names_cur_elec_bs = [names_cur_elec_bs; convertCharsToStrings(cur_name)];
                    % names_cur_elec_ind_bs = [names_cur_elec_ind_bs; elec_ind(j)];
        
                    if contains(cur_name,"'")
                        cur_num = cur_name(3:numel(cur_name));
                    else
                        cur_num = cur_name(2:numel(cur_name));
                    end
        
                    nums_cur_elec_bs = [nums_cur_elec_bs; str2double(convertCharsToStrings(cur_num))];
                end
            
            else
                %if elec_letters(i) is M but M'left also exists, then we don't want M'left to count as containing M
                if contains(cur_name,elec_letters(i)) && ~contains(cur_name,"'")
                    names_cur_elec_bs = [names_cur_elec_bs; convertCharsToStrings(cur_name)];
                    % names_cur_elec_ind_bs = [names_cur_elec_ind_bs; elec_ind(j)];
        
                    %if contains(cur_name,"'")
                    %    cur_num = cur_name(3:numel(cur_name));
                    %else
                    cur_num = cur_name(2:numel(cur_name));
                    %end
        
                    nums_cur_elec_bs = [nums_cur_elec_bs; str2double(convertCharsToStrings(cur_num))];
                end

            end
        
        end
        [cur_outer_bs, idx_cur_outer_bs] = max(nums_cur_elec_bs); % max for THIS elec letter, bs
        cur_outer_name_bs = names_cur_elec_bs(idx_cur_outer_bs);
        %outer_contacts_bs = [outer_contacts_bs; cur_outer_bs];
        %outer_contacts_ind_bs = [outer_contacts_ind_bs; idx_cur_outer_su];
    
       
        % compare largest num in channel_names_bs with largest num in elec_name
        if min([cur_outer_bs cur_outer_su]) == cur_outer_bs
            outer_contacts = [outer_contacts; cur_outer_name_bs];
        else
            outer_contacts = [outer_contacts; cur_outer_name_su];
        end
        
    end
    
    elec_name_list = [];
    for i = 1:numel(elec_name)
        cur_name = elec_name{i};
        elec_name_list = [elec_name_list; convertCharsToStrings(cur_name)];
    end
    
    % now find the indices in elec_name corresponding to the outer_contacts
    outer_contacts_ind = [];
    for i = 1:numel(outer_contacts)
        idx_contact_name = find(elec_name_list == outer_contacts(i));
        outer_contacts_ind = [outer_contacts_ind; idx_contact_name]; % USE this to access elements in elec_ind
    end
    
    % save outer_contacts and outer_contacts_ind
    final_out_dir = sprintf('/media/Data/Human_Intracranial_MAD/analysis/Cross_Freq_Coupling/%s/%s/synchrony_data/generalized_phase/extracted_phase/outer_contacts',subject_ID, reference);
    if ~exist(final_out_dir, 'dir'), mkdir(final_out_dir); end
    
    fname_names = sprintf('%s/outer_contact_names.mat',final_out_dir);
    fname_ind = sprintf('%s/outer_contacts_ind.mat',final_out_dir);
    save(fname_names,"outer_contacts")
    save(fname_ind,"outer_contacts_ind")

end % end loop through subjects

% EVERYTHING WORKS DOWN TO HERE!!








%% Iterate through channel_name and cTable
if strcmp(reference, 'Ground')
    load(sprintf('%s/%s/%s_MAD_SES%d_Setup.mat',data_base_dir,subject_ID,subject_ID,1),"elec_name","elec_ind"); % set sesnum = 1
else
    load(sprintf('%s/%s/%s_MAD_SES%d_Setup_%s.mat',data_base_dir,subject_ID,subject_ID,1,reference),"channel_name","channel_ind");
end

cTable = readtable(sprintf('/media/Data/Human_Intracranial_MAD/analysis/Location_Consistency/3D Coordinates/%s.tsv',subject_ID), "FileType","delimitedtext",'Delimiter', '\t');
channel_names_bs = cTable.Channel; % bs = brainstorm

outer_contacts_bs = []; outer_contact_ind_bs = []; % bs = brainstorm
outer_contacts_su = []; outer_contact_ind_su = []; % su = setup file

elec_letters = [];

temp_elec_nums = []; % each electrode letter at a time: store all numbers
temp_elec_nums_bs = []; % same thing for brainstorm data

% loop throuh elec_name (from setup file)
% PROBLEM: contacts are arranged in ascending order for neighbor_average
% but descending order for Ground reference
prev_name = ' '; %prev_name = elec_name{1};
for i = 1:numel(elec_name)
    cur_name = elec_name{i};
    prev_name_char = ''; % letter of the electrode
    prev_name_num = ''; % number of the contact
    if cur_name(1) ~= prev_name(1)  
        for j = 1:numel(prev_name)
            if ~isinteger(prev_name(j))
                prev_name_char = prev_name_char + prev_name(j);
            end
        end
        for j = 1:numel(prev_name)
            if isinteger(prev_name(j))
                prev_name_num = prev_name_num + prev_name(j);
            end
        end
        elec_letters = [elec_letters; convertCharsToStrings(prev_name_char)];
        temp_elec_nums = [temp_elec_nums; str2double(convertCharsToStrings(prev_name_char))]; % TODO may need to first convert prev_name_char to str
        outer_contacts_su = [outer_contacts_su; convertCharsToStrings(prev_name)];
        outer_contact_ind_su = [outer_contact_ind_su; elec_ind(i)];
    end
    % prev_name_char = ''; % letter of the electrode
    % prev_name_num = ''; % number of the contact
    % if i == numel(elec_name)
    %     for k = 1:numel(cur_name)
    %         if ~isinteger(cur_name(k))
    %             prev_name_char = prev_name_char + k;
    %         end
    %     end
    %     for k = 1:numel(prev_name)
    %         if isinteger(prev_name(k))
    %             prev_name_num = prev_name_num + k;
    %         end
    %     end
    %     elec_letters = [elec_letters; convertCharsToStrings(prev_name_char)];
    %     temp_elec_nums = [temp_elec_nums; str2double(convertCharsToStrings(prev_name_char))];
    %     outer_contacts_su = [outer_contacts_su; convertCharsToStrings(cur_name)];
    %     outer_contact_ind_su = [outer_contact_ind_su; elec_ind(i)];
    % end
    prev_name = cur_name;
end
elec_letters = unique(elec_letters);

elec_ind_bs = 1:numel(channel_names_bs);
% loop through brainstorm contacts (channel_names_bs)
prev_name = channel_names_bs{1}; % TODO may not be correct way to access
for i = 2:numel(channel_names_bs)
    cur_name = channel_names_bs{i};
    prev_name_num = ''; % number of the contact
    if cur_name(1) ~= prev_name(1)
        for j = 1:numel(prev_name)
            if isinteger(prev_name(j))
                prev_name_num = prev_name_num + j;
            end
        end
        temp_elec_nums_bs = [temp_elec_nums_bs; str2double(convertCharsToStrings(prev_name_char))]; % TODO may need to first convert prev_name_char to str
        outer_contacts_bs = [outer_contacts_bs; convertCharsToStrings(prev_name)];
        outer_contact_ind_bs = [outer_contact_ind_bs; elec_ind_bs(i-1)];
    end
    prev_name_num = ''; % number of the contact
    if i == numel(channel_names_bs)
        for k = 1:numel(prev_name)
            if isinteger(prev_name(k))
                prev_name_num = prev_name_num + k;
            end
        end
        temp_elec_nums_bs = [temp_elec_nums_bs; str2double(convertCharsToStrings(prev_name_char))];
        outer_contacts_bs = [outer_contacts_bs; convertCharsToStrings(cur_name)];
        outer_contact_ind_bs = [outer_contact_ind_bs; elec_ind_bs(i)];
    end
    prev_name = cur_name;
end

% reconcile outer contacts from bs and su
outer_contacts = []; % final version to save
outer_contacts_ind = [];

% outer_contacts_bs and outer_contacts_su should have same number of elements
for i = 1:numel(outer_contacts_bs)
    [~,min_outer_num] = min([temp_elec_nums_bs(i) temp_elec_nums(i)]);
    if min_outer_num == 1
        outer_contacts = [outer_contacts; outer_contacts_bs(i)];
        outer_contacts_ind = [outer_contacts_ind; outer_contact_ind_bs(i)];
    else
        outer_contacts = [outer_contacts; outer_contacts_su(i)];
        outer_contacts_ind = [outer_contacts_ind; outer_contact_ind_su(i)];
    end
end

% save outer_contacts and outer_contact_ind
final_out_dir = sprintf('/media/Data/Human_Intracranial_MAD/analysis/Cross_Freq_Coupling/%s/%s/synchrony_data/generalized_phase/extracted_phase/outer_contacts',subject_ID, reference);
if ~exist(final_out_dir, 'dir'), mkdir(final_out_dir); end

fname_names = sprintf('%s/outer_contact_names.mat',final_out_dir);
fname_ind = sprintf('%s/outer_contacts_ind.mat',final_out_dir);
save(fname_names,"outer_contacts")
save(fname_ind,"outer_contacts_ind")