%% leadfield check by Vincent
% 2018.4.2 modified from Pedrito's code
clc;clear;
DISK = 'F';
% Load Cortex surface and Eledtrodes Positions: name.pos
LF_PATH = strcat(DISK,':\00TychoMonkey\20170920_TychoMonkey_Vincent\res\LF\');
SV_PATH = strcat(DISK,':\00TychoMonkey\20170920_TychoMonkey_Vincent\res\LF\LF-Check\');
% load vis tool
addpath(strcat(DISK,':\00TychoMonkey\20170920_TychoMonkey_Vincent\Tools\Human-Visualization\'));
%addpath(strcat(DISK,':\00TychoMonkey\20170920_TychoMonkey_Vincent\Tools\iso2mesh\'));
load(fullfile(LF_PATH, 'Cortex_red_Su.mat')); Cortex = Cortex_red; clear Cortex_red;
% Cortex_red = reducepatch(Cortex,5e4/size(Cortex.vertices,1));
% [Cortex_red.vertices,Cortex_red.faces] = meshcheckrepair(Cortex_red.vertices,Cortex_red.faces,'meshfix');
EEG_ele  = load(fullfile(LF_PATH, 'EEG-elecs_Su.mat'));
ECoG_ele = load(fullfile(LF_PATH, 'EcoG-elecs_Su-embbed.mat'));
EEG_ele.pos  = EEG_ele.electrodes; ECoG_ele.pos = ECoG_ele.XYZ;
% load leadfield: name.K
load(fullfile(LF_PATH, 'LF\Monkey-EEG.mat'));  EEG_LF  = K; clear K;
load(fullfile(LF_PATH, 'LF\Monkey-EcoG.mat')); ECoG_LF = K; clear K;
% select the leadfield

% %% change the lead field order for EEG
% % For the i-th electrodes (you can do this for all EEG electrodes and some EcoG ones)
% names1020 = char('Fp1','Fp2','F3','F4','C3','C4','P3','P4','O1','O2','F7','F8','T3','T4','T5','T6','FZ','PZ','A1','A2');
% namesnaotaka = char('Fp1','Fp2','F7','F3','FZ','F4','F8','T3','C3','C4','T4','T5','P3','PZ','P4','T6','O1','O2');
% [~,loc] = ismember(namesnaotaka, names1020, 'rows');
% clear names1020 namesnaotaka;
% EEG_LF = EEG_LF(nonzeros(loc),:);
% save([SV_PATH, 'Monkey-EEG-18'],'EEG_LF'); save([SV_PATH, 'Monkey-ECoG-128'],'ECoG_LF'); 
%% check leadfield and topography for EEG
% For the i-th electrodes (you can do this for all EEG electrodes and some EcoG ones)

MODALITY = 'EEG';
K= EEG_LF;
[Ne, Ns3] = size(K);
for i_ele = 1:length(EEG_ele.pos)
    i_ele = 15;
    Ki = K(i_ele,:);
    Ki = reshape(Ki,[3 Ns3/3])';
    Kin = sqrt(dot(Ki,Ki,2));
    plotfield(1,1,0,'none',Cortex, Cortex.vertices, Kin, 'jet');
    colorbar;
    hold on;
    scatter3(EEG_ele.pos(i_ele,1),EEG_ele.pos(i_ele,2),EEG_ele.pos(i_ele,3),40,'k','filled');
    hold on;
    P = Cortex.vertices(1:50:end,:);
    Kj = Ki(1:50:end,:);
    quiver3(P(:,1),P(:,2), P(:,3),Kj(:,1), Kj(:,2) , Kj(:,3), 3, 'k');
    view(-97,59);
    title(strcat(MODALITY,'-',num2str(i_ele)));
    fig=gcf;
    saveas(fig,strcat(SV_PATH,MODALITY,'-',num2str(i_ele),'.jpg'));
    close;
end

%% check leadfield and topography for ECoG
MODALITY = 'ECoG';
K= ECoG_LF;
[Ne, Ns3] = size(K);
for i_ele = 1:length(ECoG_ele.pos)
    %i_ele = 12;
    Ki = K(i_ele,:);
    Ki = reshape(Ki,[3 Ns3/3])';
    Kin = sqrt(dot(Ki,Ki,2));
    plotfield(1,1,0,'none',Cortex, Cortex.vertices, Kin, 'jet');
    colorbar;
    hold on;
    P = Cortex.vertices(1:50:end,:);
    Kj = Ki(1:50:end,:);
    quiver3(P(:,1),P(:,2), P(:,3),Kj(:,1), Kj(:,2) , Kj(:,3), 5, 'k');
    view(-97,59);
    hold on;
    scatter3(ECoG_ele.pos(i_ele,1),ECoG_ele.pos(i_ele,2),ECoG_ele.pos(i_ele,3),40,'k','filled');
    title(strcat(MODALITY,'-',num2str(i_ele)));
    fig=gcf;
    saveas(fig,strcat(SV_PATH,MODALITY,'-',num2str(i_ele),'.jpg'));
    close;
end

%% check source for EEG 
% For the j-th source point (choose some extreme points to test)
MODALITY = 'EEG';
K= EEG_LF;
[Ne,Ns3] = size(K);
for j_source  =1:1000:Ns3/3
j_source = 5856; % you can decide bigger circles here
Kj = K(:,(3*(j_source-1)+1):(3*(j_source-1)+3));
Kjn = sqrt(dot(Kj,Kj,2));
plotsurf(Cortex,1,.1,1,'',[0 1 0],'r');
hold on;
scatter3(EEG_ele.pos(:,1), EEG_ele.pos(:,2), EEG_ele.pos(:,3), 100*Kjn, 'r', 'filled');
text( Cortex.vertices(j_source,1), Cortex.vertices(j_source,2), Cortex.vertices(j_source,3), 'X','color', 'k');
title(strcat(MODALITY,'-source-',num2str(j_source)));
fig=gcf;
saveas(fig,strcat(SV_PATH,MODALITY,'-source-',num2str(j_source),'.jpg'));
close;
end

%% check source for ECoG
MODALITY = 'ECoG';
K= ECoG_LF;
[Ne,Ns3] = size(K);
for j_source  =1:1000:Ns3/3
%j_source = 4556; % you can decide bigger circles here
Kj = K(:,(3*(j_source-1)+1):(3*(j_source-1)+3));
Kjn = sqrt(dot(Kj,Kj,2));
plotsurf(Cortex,1,.1,1,'',[0 1 0],'r');
hold on;
scatter3(ECoG_ele.pos(:,1), ECoG_ele.pos(:,2), ECoG_ele.pos(:,3), 10*Kjn, 'r', 'filled');
text( Cortex.vertices(j_source,1), Cortex.vertices(j_source,2), Cortex.vertices(j_source,3), 'X','color', 'k');
title(strcat(MODALITY,'-source-',num2str(j_source)));
fig=gcf;
saveas(fig,strcat(SV_PATH,MODALITY,'-source-',num2str(j_source),'.jpg'));
close;
end

%%%%%%%%%%%%%%%%%%%%%%
% %% old code 
% plotsurf(Inskull,1,.1,1,'',[0 1 0],'r')
% hold on; scatter3(pos2(1:10:end,1),pos2(1:10:end,2),pos2(1:10:end,3),30,'b','filled');
% plotsurf(Silicone,0,.1,1,'',[1 1 0],'r')
% hold on; scatter3(pos1(1:10:end,1),pos1(1:10:end,2),pos1(1:10:end,3),30,'r','filled');
% face= Scalp.faces; node = Scalp.vertices
% trisurf(Scalp.faces(:,1:3),Scalp.vertices(:,1),Scalp.vertices(:,2),Scalp.vertices(:,3),'color','r')
% hold on; scatter3(electrodes(:,1),electrodes(:,2),electrodes(:,3),30,'b','filled');
%%%%%%%%%%%%%%%%%%%
% %% old version 
% clear; clc; close all;
% restoredefaultpath;
% DISK = 'G';
% proj_path = strcat(DISK,':\00TychoMonkey\20170920_TychoMonkey_Vincent\');
% 
% load([proj_path,  'res\leadfieldSurface\ecog_eeg_lead_fields\Monkey-EEG.mat']); K_EEG  = K; clear K;
% load([proj_path,  'res\leadfieldSurface\ecog_eeg_lead_fields\Monkey-EcoG.mat']); K_EcoG = K; clear K;
% load([proj_path,  'res\leadfieldSurface\ecog_eeg_lead_fields\Su_cortex5232_normals.mat']);
% load([proj_path,  'res\leadfieldSurface\ecog_eeg_lead_fields\Su_cortex5232_trian.mat']);
% load([proj_path,  'res\leadfieldSurface\ecog_eeg_lead_fields\Su_cortex5232_coo.mat']);
% 
% cortex_test.faces = cortex5000_trian;
% %cortex_test.FaceNormals = cortex_normals;
% cortex_test.vertices = cortex5000_coo;
% 
% KK = kron(eye(5232),[1,1,1]');
% TR=K_EEG*KK;
% 
% cortex_opt.texture = ones(length(cortex_test.vertices),1);
% spm_eeg_render(cortex_test,cortex_opt);
% cortex_opt.texture = TR(70,:)';
% %cortex_opt.texture = ecog_LF(89,:)';
% spm_eeg_render(cortex_test,cortex_opt);

