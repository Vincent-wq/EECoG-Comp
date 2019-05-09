function [data_EEG, data_ECoG] = ft_format_MonkeyEEG_ECoG_Data(dat_path, dat_sv_path)
%% reading the original and down-sampled Simutaneuous EEG/ECoG experiment data
% EEG and ECoG data
% By Qing Wang (Vincent) 11th Oct. 2017
% dat_path = 'G:\00TychoMonkey\20170918_originalTychoMonkey_Su\EEG4096_ECoG1000';
%dat_path = 'F:\00TychoMonkey\20170918_originalTychoMonkey_Su\EEG1000_ECoG1000';
addpath(genpath(dat_path));
% Constants data type
dat_type = dat_path(end-12:end-9);
% loading data and format to cell;
if strcmp(dat_type, '4096')
    disp('loading raw data with fs_EEG = 4096  Hz and fs_ECoG = 1000  Hz...')
    Fs_EEG = 4096;
    % EEG data
    % awake
    awake01_EEG_data = load('EEG01.mat'); awake02_EEG_data = load('EEG02.mat');
    % anetheisa
    anesthesia01_EEG_data = load('EEG01_anesthesia.mat'); anesthesia02_EEG_data = load('EEG02_anesthesia.mat');
    anesthesia03_EEG_data = load('EEG03_anesthesia.mat'); anesthesia04_EEG_data = load('EEG04_anesthesia.mat');
    anesthesia05_EEG_data = load('EEG05_anesthesia.mat'); 
    EEG_raw_all = {awake01_EEG_data.EEG_rest_01, awake02_EEG_data.EEG_rest_02, ...
        anesthesia01_EEG_data.EEG_anes_01, anesthesia02_EEG_data.EEG_anes_02,  ...
        anesthesia03_EEG_data.EEG_anes_03, anesthesia04_EEG_data.EEG_anes_04,  ...
        anesthesia05_EEG_data.EEG_anes_05};
    % ECoG data
    % awake
    awake01_ECoG_data = load('ECoG01.mat'); awake02_ECoG_data = load('ECoG02.mat');
    % anetheisa
    anesthesia01_ECoG_data = load('ECoG01_anesthesia.mat'); anesthesia02_ECoG_data = load('ECoG02_anesthesia.mat');
    anesthesia03_ECoG_data = load('ECoG03_anesthesia.mat'); anesthesia04_ECoG_data = load('ECoG04_anesthesia.mat');
    anesthesia05_ECoG_data = load('ECoG05_anesthesia.mat'); 
    ECoG_raw_all = {awake01_ECoG_data.WaveData, awake02_ECoG_data.WaveData, ...
        anesthesia01_ECoG_data.WaveData, anesthesia02_ECoG_data.WaveData,   ...
        anesthesia03_ECoG_data.WaveData, anesthesia04_ECoG_data.WaveData,   ...
        anesthesia05_ECoG_data.WaveData}; 
elseif strcmp(dat_type, '1000')
    disp('loading down sampled data with fs_EEG = 1000Hz and fs_ECoG =  1000Hz...')
    Fs_EEG = 1000;
    % EEG data
    % awake
    awake01_EEG_data = load('EEG01.mat'); awake02_EEG_data = load('EEG02.mat');
    % anetheisa
    anesthesia01_EEG_data = load('EEG01_anesthesia.mat'); anesthesia02_EEG_data = load('EEG02_anesthesia.mat');
    anesthesia03_EEG_data = load('EEG03_anesthesia.mat'); anesthesia04_EEG_data = load('EEG04_anesthesia.mat');
    anesthesia05_EEG_data = load('EEG05_anesthesia.mat'); 
    EEG_raw_all = {awake01_EEG_data.EEG2, awake02_EEG_data.EEG2,...
        anesthesia01_EEG_data.EEG2, anesthesia02_EEG_data.EEG2,...
        anesthesia04_EEG_data.EEG2, anesthesia05_EEG_data.EEG2,...
        anesthesia01_EEG_data.EEG2};
    % ECoG data
    % awake
    awake01_ECoG_data = load('ECoG01.mat'); awake02_ECoG_data = load('ECoG02.mat');
    % anetheisa
    anesthesia01_ECoG_data = load('ECoG01_anesthesia.mat'); anesthesia02_ECoG_data = load('ECoG02_anesthesia.mat');
    anesthesia03_ECoG_data = load('ECoG03_anesthesia.mat'); anesthesia04_ECoG_data = load('ECoG04_anesthesia.mat');
    anesthesia05_ECoG_data = load('ECoG05_anesthesia.mat'); 
    ECoG_raw_all = {awake01_ECoG_data.WaveData, awake02_ECoG_data.WaveData,...
        anesthesia01_ECoG_data.WaveData, anesthesia02_ECoG_data.WaveData,...
        anesthesia03_ECoG_data.WaveData, anesthesia04_ECoG_data.WaveData,...
        anesthesia05_ECoG_data.WaveData};
else
    disp('Error£º unrecognized data ...')
end
% clear unused data
clear *_EEG_data; clear *_ECoG_data;  clear EEG0*; clear WaveData;

%% fromating FieldTrip trial Label
N_trial_EEG = length(EEG_raw_all); N_trial_ECoG = length(ECoG_raw_all); 
TRIAL_NAME = {'awake1', 'awake2', 'anesthesia1', 'anesthesia2', ...
    'anesthesia3', 'anesthesia4', 'anesthesia5' };

% Sampling frequency
Fs_ECoG =1000;

% label Names
%EEG
EEG_label = {'Fp1','Fp2','F7','F3','Fz','F4','F8','T3','C3','C4','T4','T5','P3','Pz','P4','T6','O1','O2','Trigger'};

%ECoG
ECoG_label = [];
for ilabel = 1:128
    ECoG_label{ilabel} = strcat('Ch_',num2str(ilabel));
end
ECoG_label{129} = 'Trigger';

% display the information of synchronization
info_Flag = 1;

% Formating
for i_eeg = 1 : N_trial_EEG 
    eeg_tmp = {};
    % sampling frequency
    eeg_tmp.fsample = Fs_EEG;
    % channel label
    eeg_tmp.label = EEG_label;
    % trial
    eeg_tmp.trial= {double(ft_synchronize_trigger_for_trial(EEG_raw_all{i_eeg}, Fs_EEG, info_Flag))};
    % trl information
    trl_EEG = [];
    trl_EEG(1) = 1; trl_EEG(2) = length(eeg_tmp.trial{1}); trl_EEG(3) = 0;
    % trial time 
    eeg_tmp.time = {(1:size(eeg_tmp.trial{1},2))/eeg_tmp.fsample};
    % formatting
    cfg_EEG = [];
    cfg_EEG.trl = trl_EEG;
    eeg_tmp = ft_redefinetrial(cfg_EEG, eeg_tmp);
    % save data to one variable
    eval(strcat('data_EEG.', TRIAL_NAME{i_eeg}, '=','eeg_tmp;'))
end

for i_ecog = 1 : N_trial_ECoG
    ecog_tmp = {};
    % sampling frequency
    ecog_tmp.fsample = Fs_ECoG;
    % channel label
    ecog_tmp.label = ECoG_label;
    % trial
    ecog_tmp.trial = {double(ft_synchronize_trigger_for_trial(ECoG_raw_all{i_ecog}, Fs_ECoG, info_Flag))};
    % trl information
    trl_ECoG = [];
    trl_ECoG(1) = 1; trl_ECoG(2) = length(ecog_tmp.trial{1}); trl_ECoG(3) = 0;
    % trial time 
    ecog_tmp.time = {(1:size(ecog_tmp.trial{1},2))/ecog_tmp.fsample};
    % formatting
    cfg_ECoG = [];
    cfg_ECoG.trl = trl_ECoG;
    ecog_tmp = ft_redefinetrial(cfg_ECoG, ecog_tmp);
    % save data to one variable
    eval(strcat('data_ECoG.', TRIAL_NAME{i_ecog}, '=','ecog_tmp;'))
end

clear the tmp data:
clear cfg*;  clear trl*; clear N*; clear *raw_all;
% save data 
if ~isempty(dat_sv_path)
  sv_path = dat_sv_path;
  if strcmp(dat_type, '1000')
      save(strcat(sv_path, 'ft_data_sync_EEG1000'), 'data_EEG');
      save(strcat(sv_path, 'ft_data_sync_ECoG1000'), 'data_ECoG', '-v7.3');
  elseif strcmp(dat_type, '4096')
      save(strcat(sv_path, 'ft_data_sync_EEG4096r'), 'data_EEG');
      save(strcat(sv_path, 'ft_data_sync_ECoG1000r'), 'data_ECoG', '-v7.3');
  end
end
data_EEG; data_ECoG;


