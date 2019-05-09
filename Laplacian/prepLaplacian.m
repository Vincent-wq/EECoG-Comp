clc; clear; close all;
P= cd; FolderLength = length('1-Laplacian');
addpath(genpath([P(1:end-FolderLength) 'tools\chronux_2_11']));
addpath(genpath([P(1:end-FolderLength) 'tools\RobustSpectra']));
addpath([P(1:end-FolderLength) '\tools\eeglab14_1_2b\']);
addpath(genpath([P(1:end-FolderLength),'\tools\CSDtoolbox\']));

addpath(genpath([P(1:end-FolderLength),'\dataWD\']));
%% Calculate EEG Laplacian 
eegNames={'eegAwake-500','eegAwakeAr-500','eegAne-500','eegAneAr-500','Cuba-200'};
neegDat = length(eegNames); 
% prepare montage and conf.s for EEG
% E = textread('ChanName15.txt','%s');                                      
% M = ExtractMontage('10-5-System_Mastoids_EGI129.csd',E);  
% [ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;
% for i = 1:neegDat
%     EEG = pop_loadset( 'filename',[eegNames{i} '.set']);
%     [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
%     [ALLEEG, EEG, CURRENTSET] = calLaplacianEEGLAB(EEG, ALLEEG, M, 0);
%     NAME = [eegNames{i}, '-Lap'];
%     pop_saveset(EEG, 'filename',[NAME '.set']);
% end
% eeglab redraw

%% Calculate ECoG Laplacian 
ecogNames={'ecogAwake-500','ecogAwakeAr-500','ecogAne-500','ecogAneAr-500'};
necogDat = length(ecogNames);
% calculate Laplacian on the cortex
ecogEleF = 'ECoG_ele128'; cortexF = 'Cortex_red_10007'; geoFName = 'geoLp';
%[geoFile]=ecogGeo4Laplacian(ecogEleF,cortexF,geoFName);
load(geoFName);
% Mats = zeros(nsL,nChan,necogDat);
Mats=load('H','Mats'); Mats=Mats.Mats;
% read data
[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;
for i = 1:necogDat
    EEG = pop_loadset( 'filename',[ecogNames{i} '.set']);
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
    Vecog = EEG.data; Vecog=Vecog(:,1:end-1);
%     [lambda, LAMBDA, GCV] = LambdaGCV(S, L, Vecog, 1000, -1, 5, [ecogNames{i} '-lambda' ]);
%     [H] = calSourLap(S, L, lambda);
    H = Mats(:,:,i);
    LecogS = H*Vecog;                                % Source Laplacian
    Lecog = LecogS(inds(j),:);
    [ALLEEG, EEG, CURRENTSET] = ecogLaplacianFormat(EEG, ALLEEG, Lecog);
    NAME = [ecogNames{i}, '-Lap'];
    pop_saveset(EEG, 'filename',[NAME '.set']);
end
save('H','Mats','ecogNames');
eeglab redraw
% estimate lambda for smoothness
% reconstruct the signal on the source