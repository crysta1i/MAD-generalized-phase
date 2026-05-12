% Run this function at the start of all scripts 

function [datapath, data_base_dir, tw_out_dir, out_dir, num_sessions, Fs, elec_letters] = tw_setup(subject_ID, reference)

datapath = '/media/Data/Human_Intracranial_MAD/';
data_base_dir = sprintf('%s1_formatted',datapath);
tw_out_dir = sprintf('/media/Data/Human_Intracranial_MAD/analysis/TravWaves');
out_dir = sprintf('/media/Data/Human_Intracranial_MAD/analysis/Cross_Freq_Coupling/%s/%s/synchrony_data/generalized_phase',subject_ID, reference); 

addpath('/media/Data/Human_Intracranial_MAD/analysis/PAC_code/CF_Coupling/generalized phase/circstat-matlab')
addpath('/media/Data/Human_Intracranial_MAD/analysis/PAC_code/CF_Coupling/generalized phase/generalized-phase')
addpath('/media/Data/Human_Intracranial_MAD/analysis/PAC_code/CF_Coupling/generalized phase/generalized-phase/analysis')
addpath('/media/Data/Human_Intracranial_MAD/_toolbox/spectral-analysis')
addpath('/media/Data/Human_Intracranial_MAD/_toolbox')

addpath('/media/Data/Human_Intracranial_MAD/analysis/TravWaves/Code')
addpath('/media/Data/Human_Intracranial_MAD/analysis/TravWaves/Code/generalized-phase')

switch subject_ID
    case 'EMU001'
        num_sessions = 3;
    case 'EMU024'
        num_sessions = 3;
    case 'EMU025'
        num_sessions = 2;
    case 'EMU030'
        num_sessions = 2;
    case 'EMU037'
        num_sessions = 4;
    case 'EMU038'
        num_sessions = 1;
    case 'EMU039'
        num_sessions = 4;
    case 'EMU041'
        num_sessions = 9;
    case 'EMU047'
        num_sessions = 1;
    case 'EMU051'
        num_sessions = 1;
end

if strcmp(subject_ID,'EMU001')
    Fs = 1000;
else
    Fs = 2048;
end

if strcmp(reference,'Ground')
    load(sprintf('%s/%s/%s_MAD_SES%d_Setup.mat',data_base_dir,subject_ID,subject_ID,1),"elec_name");
else
    load(sprintf('%s/%s/%s_MAD_SES%d_Setup_%s.mat',data_base_dir,subject_ID,subject_ID,1,reference),"channel_name");
    elec_name = channel_name;
end

% get elec letters
elec_letters = [];
for i = 1:numel(elec_name)
    cur_name = elec_name{i};
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

% subject_IDs = {'EMU001','EMU024','EMU025','EMU030','EMU037','EMU038','EMU039','EMU041','EMU047','EMU051'};
% subject_IDs_str = ["EMU001","EMU024","EMU025","EMU030","EMU037","EMU038","EMU039","EMU041","EMU047","EMU051"];
% 
% num_sess_dict = dictionary(subject_IDs_str, [3, 3, 2, 2, 4, 1, 4, 9, 1, 1]);
% Fs_dict = dictionary(subject_IDs_str, [1000, 2048, 2048, 2048, 2048, 2048, 2048, 2048, 2048, 2048]);

end