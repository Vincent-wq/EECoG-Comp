% This is the main exploratory code for the simultaneous EEG and ECoG 
% neurotycho dataset (http://neurotycho.org/); 
% started from 2017.9.29 by Vincent 
%% Preparing working space
clear; clc; close all;
restoredefaultpath;
DISK = 'F';
proj_path = strcat(DISK,':\00TychoMonkey\20170920_TychoMonkey_Vincent\');
path_sv  = [proj_path, 'res\'];
path_LF = [proj_path, 'res\leadfieldSurface\'];
addpath(proj_path); addpath(path_sv);
addpath([proj_path,'Bosch']);
addpath([proj_path,'PedritoResults\5-Leadfields']);

addpath([proj_path,'Tools\fieldtrip-20170925']);
addpath(genpath([proj_path,'Tools\chronux\chronux_2_12']));
addpath([proj_path,'Tools\TARA_software']);
addpath([proj_path,'Tools\GraphVar-1.03']);
addpath([proj_path,'Tools\HERMES Toolbox 2017-01-31']);
addpath([proj_path,'Tools\Brainstorm3']);
cd(proj_path);
%% 0. Prepare all the data used in this script.
% data preprocessing 
% ft_prepare_ALL_data(proj_path);
%% 1. Load Preprocessed data
% findings: the artifact apears in all the EEG only apear on one channel (93) of ECoG
% load best EEG awake / ane data
% load('ft_dn_notch_hp_rsmp_EEG_awake1.mat'); eeg_awake1 = ft_dn_notch_hp_rsmp; clear ft_dn_notch_hp_rsmp;
% load('ft_dn_notch_hp_rsmp_EEG_anesthesia5.mat'); eeg_ane5 = ft_dn_notch_hp_rsmp; clear ft_dn_notch_hp_rsmp;
% % load best ECoG awake / ane data
% load('ft_dn_notch_hp_rsmp_ECoG_awake1.mat'); ecog_awake1 = ft_dn_notch_hp_rsmp; clear ft_dn_notch_hp_rsmp;
% load('ft_dn_notch_hp_rsmp_ECoG_anesthesia5.mat'); ecog_ane5 = ft_dn_notch_hp_rsmp; clear ft_dn_notch_hp_rsmp;
% 
% load best EEG awake / ane data
load('ft_dn_notch_hp_rsmp_EEG_awake1.mat'); eeg_awake1 = ft_dn_notch_hp_rsmp; clear ft_dn_notch_hp_rsmp;
load('ft_dn_notch_hp_rsmp_EEG_anesthesia5.mat'); eeg_ane5 = ft_dn_notch_hp_rsmp; clear ft_dn_notch_hp_rsmp;
% load best ECoG awake / ane data
load('ft_dn_notch_hp_rsmp_ECoG_awake1.mat'); ecog_awake1 = ft_dn_notch_hp_rsmp; clear ft_dn_notch_hp_rsmp;
load('ft_dn_notch_hp_rsmp_ECoG_anesthesia5.mat'); ecog_ane5 = ft_dn_notch_hp_rsmp; clear ft_dn_notch_hp_rsmp;
% 
% % load original awake /ane data
% 
% load('F:\00TychoMonkey\20170920_TychoMonkey_Vincent\res\ft_sync_EEG4096.mat');
% load('F:\00TychoMonkey\20170920_TychoMonkey_Vincent\res\ft_sync_ECoG1000.mat')
% eeg_awake1_o = data_EEG.awake1; eeg_ane5_o = data_EEG.anesthesia5; clear data_EEG;
% ecog_awake1_o = data_ECoG.awake1; ecog_ane5_o = data_ECoG.anesthesia5; clear data_ECoG;
% 
% % visual check
% cfg = []; cfg.viewmode = 'vertical'; cfg.continuous = 'yes'; cfg.ylim = 'maxmin';
% cfg.blocksize = 90; cfg.verticalpadding = 'auto';
% ft_databrowser(cfg, eeg_ane5); ft_databrowser(cfg, eeg_ane5_o);
% ft_databrowser(cfg, ecog_ane5);   ft_databrowser(cfg, ecog_ane5_o);


%% load old lead field and check
% load([path_LF, 'Cortex-mid_Su_reduced']);  cortex_mid = nP;   clear nP;
% load([path_LF, 'EEG_LF']); EEG_LF = eeg_LF; clear eeg_LF;
% load([path_LF, 'ECoG_LF']); ECoG_LF = ecog_LF; clear ecog_LF;
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
% 
% for i = 1:74
%     cortex_opt.texture = TR(i,:)';
%     spm_eeg_render(cortex_test,cortex_opt);
%     title(num2str(i))
%     caxis([0 0.5])
%     saveas(gcf,[path_sv,num2str(i)],'jpg');
%     close all;
% end
% cortex_opt_eeg.texture = sqrt(sum(eeg_LF.^2));
% spm_eeg_render(cortex_mid,cortex_opt_eeg);
% hist(cortex_opt_eeg.texture,100)
% 
% cortex_opt_ecog.texture = sqrt(sum(ecog_LF.^2));
% spm_eeg_render(cortex_mid,cortex_opt_ecog);


%% the spectrgram for ft data (4 channel).
ALL1 = [0.3 100]; ALL2    = [0.3 45]; DELTA    = [0.3 4]; THETA = [4 8]; 
ALPHA = [8, 13] ; LOWBETA = [14 20] ; HIGHBETA = [20 30]; GAMMA = [31 45]; 

ft_mtspecgramc(eeg_awake1,  1000, [4 7], ALL2, [2, 1], '', 'on');
ft_mtspecgramc(eeg_ane5,    1000, [4 7], ALL2, [2, 1], '', 'on');
ft_mtspecgramc(ecog_awake1, 1000, [4 7], ALL2, [2, 1], '', 'on');
ft_mtspecgramc(ecog_ane5,   1000, [4 7], ALL2, [2, 1], '', 'on');

%% coherence
data1 =  eeg_awake1.trial{1}; data2 =  eeg_ane5.trial{1}; 
data3 = ecog_awake1.trial{1}; data4 = ecog_ane5.trial{1}; 
win = [2 2]; 
params.Fs = 1000; params.tapers = [3 5]; params.fpass = ALL2; 

[Sc,Cmat,Ctot,Cvec,Cent,f]   = CrossSpecMatc(data1',5, params);

[Sc1,Cmat1,Ctot,Cvec,Cent,f] = CrossSpecMatc(data3',5,params);
i = 5;
h=figure('visible','on','PaperSize',[21 3]);
imagesc(abs(squeeze(Cmat1(i,:,:))));
colorbar()
caxis([0 1])
colormap jet
title(strcat('f =', num2str(f(i)), ' Hz'))

% ploting the cross-spectrum
[Sc,Cmat,Ctot,Cvec,Cent,f] = CrossSpecMatc(data1',5, params);
for i = 1:length(f)
h=figure('visible','off');
imagesc(log(abs(squeeze(Sc(i,:,:)))));
colorbar()
caxis([0 10])
colormap jet
xticks(linspace(1,19,19))
xticklabels(eeg_awake1.label)
yticks(linspace(1,19,19))
yticklabels(eeg_awake1.label)
title(strcat('f =',num2str(f(i)),' Hz'))
saveas(h, strcat('fig\cs\eeg_awake1_cs_',num2str(f(i)),'.jpg'))
imagesc(abs(squeeze(Cmat(i,:,:))));
colorbar()
caxis([0 1])
colormap jet
xticks(linspace(1,19,19))
xticklabels(eeg_awake1.label)
yticks(linspace(1,19,19))
yticklabels(eeg_awake1.label)
title(strcat('f =',num2str(f(i)),' Hz'))
saveas(h, strcat('fig\cs\eeg_awake1_coh_',num2str(f(i)),'.jpg'))
end

[Sc,Cmat,Ctot,Cvec,Cent,f] = CrossSpecMatc(data2',5,params);
for i = 1:length(f)
h=figure('visible','off');
imagesc(log(abs(squeeze(Sc(i,:,:)))));
colorbar()
caxis([0 10])
colormap jet
xticks(linspace(1,19,19))
xticklabels(eeg_awake1.label)
yticks(linspace(1,19,19))
yticklabels(eeg_awake1.label)
title(strcat('f =',num2str(f(i)),' Hz'))
saveas(h, strcat('fig\cs\eeg_ane5_cs_',num2str(f(i)),'.jpg'))
imagesc(abs(squeeze(Cmat(i,:,:))));
colorbar()
caxis([0 1])
colormap jet
xticks(linspace(1,19,19))
xticklabels(eeg_awake1.label)
yticks(linspace(1,19,19))
yticklabels(eeg_awake1.label)
title(strcat('f =',num2str(f(i)),' Hz'))
saveas(h, strcat('fig\cs\eeg_ane5_coh_',num2str(f(i)),'.jpg'))
end

[Sc,Cmat,Ctot,Cvec,Cent,f] = CrossSpecMatc(data3',5, params);
for i = 1:length(f)
h=figure('visible','off');
imagesc(log(abs(squeeze(Sc(i,:,:)))));
colorbar()
caxis([0 10])
colormap jet
title(strcat('f =',num2str(f(i)),' Hz'))
saveas(h, strcat('fig\cs\ecog_awake1_cs_',num2str(f(i)),'.jpg'))
imagesc(abs(squeeze(Cmat(i,:,:))));
colorbar()
caxis([0 1])
colormap jet
title(strcat('f =',num2str(f(i)),' Hz'))
saveas(h, strcat('fig\cs\ecog_awake1_coh_',num2str(f(i)),'.jpg'))
end

[Sc,Cmat,Ctot,Cvec,Cent,f] = CrossSpecMatc(data4',5,params);
for i = 1:length(f)
h=figure('visible','off');
imagesc(log(abs(squeeze(Sc(i,:,:)))));
colorbar()
caxis([0 10])
colormap jet
title(strcat('f =',num2str(f(i)),' Hz'))
saveas(h, strcat('fig\cs\ecog_ane5_cs_',num2str(f(i)),'.jpg'))
imagesc(abs(squeeze(Cmat(i,:,:))));
colorbar()
caxis([0 1])
colormap jet
title(strcat('f =',num2str(f(i)),' Hz'))
saveas(h, strcat('fig\cs\ecog_ane5_coh_',num2str(f(i)),'.jpg'))
end
% ploting the coherency
params.err = [2 0.01];
[Cmn,Phimn,Smn,Smm,f,ConfC,PhiStd,Cerr] = coherencyc_unequal_length_trials([data1';data2'],[5 5], params, [1 size(data1,2);size(data1,2)+1 size(data1,2)+size(data2,2) ]);
for i = 1:length(f)
h=figure('visible','off');
imagesc(abs(squeeze(Sc(i,:,:))));
colorbar()
saveas(h, strcat('fig\cs\ecog_ane5_cs_',num2str(f(i)),'.jpg'))
imagesc(abs(squeeze(Cmat(i,:,:))));
saveas(h, strcat('fig\cs\ecog_ane5_coh_',num2str(f(i)),'.jpg'))
end


%% checking the parame chronux multitapper function
% ft_mtspecgramc(eeg_awake1, 1000, [3 5], [0 45], [2, 2], 'fig\eeg_awake1_35.jpg', 'off')
% ft_mtspecgramc(eeg_awake1, 1000, [4 7], [0 45], [2, 2], 'fig\eeg_awake1_47.jpg', 'off')
% ft_mtspecgramc(eeg_awake1, 1000, [5 9], [0 45], [2, 2], 'fig\eeg_awake1_59.jpg', 'off')
% 
% ft_mtspecgramc(eeg_ane5, 1000, [3 5], [0 45], [2, 2], 'fig\eeg_ane5_35.jpg', 'off')
% ft_mtspecgramc(eeg_ane5, 1000, [4 7], [0 45], [2, 2], 'fig\eeg_ane5_47.jpg', 'off')
% ft_mtspecgramc(eeg_ane5, 1000, [5 9], [0 45], [2, 2], 'fig\eeg_ane5_59.jpg', 'off')
% 
% ft_mtspecgramc(ecog_awake1, 1000, [3 5], [0 45], [2, 2], 'fig\ecog_awake1_35.jpg', 'off')
% ft_mtspecgramc(ecog_awake1, 1000, [4 7], [0 45], [2, 2], 'fig\ecog_awake1_47.jpg', 'off')
% ft_mtspecgramc(ecog_awake1, 1000, [5 9], [0 45], [2, 2], 'fig\ecog_awake1_59.jpg', 'off')
% 
% ft_mtspecgramc(ecog_ane5, 1000, [3 5], [0 45], [2, 2], 'fig\ecog_ane5_35.jpg', 'off')
% ft_mtspecgramc(ecog_ane5, 1000, [4 7], [0 45], [2, 2], 'fig\ecog_ane5_47.jpg', 'off')
% ft_mtspecgramc(ecog_ane5, 1000, [5 9], [0 45], [2, 2], 'fig\ecog_ane5_59.jpg', 'off')


%  'Fp1'    'Fp2'    'F3 '    'F4 '    'C3 '    'C4 '    'P3 '    'P4 '
%  'O1 '    'O2 '    'F7 '    'F8 '    'T3 '    'T4 ' 'T5 '    'T6 '    'FZ '    'PZ '

% new_label = {'Fp1','Fp2','F7','F3','Fz','F4','F8','T3','C3','C4','T4','T5','P3','Pz','P4','T6','O1','O2','Trigger'};

% exploring channel data with pburg
channelName = 'O2';
Fspan=linspace(0.01, 120, 1000);
FS=1000;
ORDER=256;
[sigy,sigx] =  ft_pburg(old_detr.awake_detr.eeg, channelName, ORDER, FS, Fspan);
[sigy,sigx] =  ft_pburg(old_detr.eeg_awake_detr_3_40, channelName, ORDER, FS, Fspan);
[sigy,sigx] =  ft_pburg(new_detr.dat_awake,    channelName, ORDER, FS, Fspan);

[sigy,sigx] =  ft_pburg(old_detr.awake_detr.eeg, channelName, ORDER, FS, Fspan);
[sigy,sigx] =  ft_pburg(new_detr.dat_awake,  channelName, ORDER, FS, Fspan);

%% using chronux 

%% Building network

%% Network Analytics

%% old codes



%% 