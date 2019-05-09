function [] = ft_prepare_ALL_data(proj_path)
% Prepare all the data used in the analysis
% by Vincent

DISK = 'F';
proj_path = strcat(DISK,':\00TychoMonkey\20170920_TychoMonkey_Vincent\');
cd(proj_path);
path_sv  = [proj_path, 'res\'];
TARApath = [proj_path,'Tools\TARA_software'];

%% 1. synchronize, formate and save the original data
% [data_EEG1000, data_ECoG1000] = ft_format_MonkeyEEG_ECoG_Data(...
%     'F:\00TychoMonkey\20170918_originalTychoMonkey_Su\EEG1000_ECoG1000', path_sv);
% [data_EEG4096, ~] = ft_format_MonkeyEEG_ECoG_Data(...
%     'F:\00TychoMonkey\20170918_originalTychoMonkey_Su\EEG4096_ECoG1000', path_sv);
% % check the data quality: synchronize
% % cfg = []; cfg.viewmode='vertical';
% % ft_databrowser(cfg, data_ECoG.anesthesia5);
% % EEG:  awake1 ok, awake2: no 39592 to end ; ane1 ok; ane2 ok; ane3 ok; ane4: no; ane5: ok;
% % ECoG: awake1 ok, awake2: ok; ane1 ok; ane2 ok; ane3 ok; ane4: ok; ane5: ok;
% % manual fix
% % EEG awake 2
% aw2 = data_EEG.awake2;
% aw2.trial = {aw2.trial{1}(:,39592:end)};
% trl_EEG = [];
% trl_EEG(1) = 1; trl_EEG(2) = length(aw2.trial{1}); trl_EEG(3) = 0;
% aw2.time = {(1:size(aw2.trial{1},2))/aw2.fsample};
% cfg_EEG = []; cfg_EEG.trl = trl_EEG;
% aw2 = ft_redefinetrial(cfg_EEG, aw2);
% data_EEG.awake2=aw2;
% % EEG awake 4
% an4 = data_EEG.anesthesia4;
% an4.trial = {an4.trial{1}(:,49874:end)};
% trl_EEG = [];
% trl_EEG(1) = 1; trl_EEG(2) = length(an4.trial{1}); trl_EEG(3) = 0;
% an4.time = {(1:size(an4.trial{1},2))/an4.fsample};
% cfg_EEG = []; cfg_EEG.trl = trl_EEG;
% an4 = ft_redefinetrial(cfg_EEG, an4);
% data_EEG.anesthesia4=an4;
% 
% save(strcat(path_sv, 'ft_sync_EEG4096'), 'data_EEG');
% save(strcat(path_sv, 'ft_sync_ECoG1000'), 'data_ECoG', '-v7.3');
% cfg = []; cfg.viewmode='vertical';
% ft_databrowser(cfg, data_EEG.anesthesia4);
% final check sync
% EEG:  awake1 ok, awake2: ok; ane1 ok; ane2 ok; ane3 ok; ane4: ok; ane5: ok;
% ECoG: awake1 ok, awake2: ok; ane1 ok; ane2 ok; ane3 ok; ane4: ok; ane5: ok;
%% 2. TARA denosing and reformat;
% take only awake2 and anesthesia5 as example to denoise
% 1. de-median; 2. TARA; 3. Notch, 50+-1Hz; 4. High pass >0.1
load([path_sv, 'ft_sync_EEG4096']);
load([path_sv, 'ft_sync_ECoG1000']);

highpassB = 0.3; notchB = 50; Nit = 50;

% EEG 
ft_preproc_sv(data_ECG.awake1, highpassB, notchB, Nit, 'EEG_awake1', path_sv, TARApath);
ft_preproc_sv(data_EEG.awake2, highpassB, notchB, Nit, 'EEG_awake2', path_sv, TARApath);
ft_preproc_sv(data_EEG.anesthesia1, highpassB, notchB, Nit, 'EEG_anesthesia1', path_sv, TARApath);
ft_preproc_sv(data_EEG.anesthesia2, highpassB, notchB, Nit, 'EEG_anesthesia2', path_sv, TARApath);
ft_preproc_sv(data_EEG.anesthesia3, highpassB, notchB, Nit, 'EEG_anesthesia3', path_sv, TARApath);
ft_preproc_sv(data_EEG.anesthesia4, highpassB, notchB, Nit, 'EEG_anesthesia4', path_sv, TARApath);
ft_preproc_sv(data_EEG.anesthesia5, highpassB, notchB, Nit, 'EEG_anesthesia5', path_sv, TARApath);

% ECoG
ft_preproc_sv(data_ECoG.awake1, highpassB, notchB, Nit, 'ECoG_awake1', path_sv, TARApath);
ft_preproc_sv(data_ECoG.awake2, highpassB, notchB, Nit, 'ECoG_awake2', path_sv, TARApath);
ft_preproc_sv(data_ECoG.anesthesia1, highpassB, notchB, Nit, 'ECoG_anesthesia1', path_sv, TARApath);
ft_preproc_sv(data_ECoG.anesthesia2, highpassB, notchB, Nit, 'ECoG_anesthesia2', path_sv, TARApath);
ft_preproc_sv(data_ECoG.anesthesia3, highpassB, notchB, Nit, 'ECoG_anesthesia3', path_sv, TARApath);
ft_preproc_sv(data_ECoG.anesthesia4, highpassB, notchB, Nit, 'ECoG_anesthesia4', path_sv, TARApath);
ft_preproc_sv(data_ECoG.anesthesia5, highpassB, notchB, Nit, 'ECoG_anesthesia5', path_sv, TARApath);

%% test code
cfg = [];
cfg.viewmode = 'vertical';
ft_databrowser(cfg, ft_dn_notch_hp_rsmp);

% prepare 
DATA = data_EEG.awake2;
Fspan=linspace(0.01, 120, 1000);
ORDER=128;
trimN=3; fs = DATA.fsample; highpassB = 0.1; notchB = 50; Nit = 50; res = 0.1;
F_Notch = designfilt('bandstopiir','FilterOrder',10, ...
         'HalfPowerFrequency1',notchB-1,'HalfPowerFrequency2',notchB+1, ...
         'SampleRate',fs);
f_hp = highpassB/(fs/2);
dHP  = fdesign.highpass('Fst,Fp,Ast,Ap', f_hp*1, f_hp*1.5, 60, 1);
F_HP = design(dHP,'butter');
% 
y = DATA.trial{1}(2,:);
y_md = y -median(y);
[x1, x2, f, cost] = tara2_L1(y_md, 2, 0.03, 0.05, 1.4, 100, 20);
f = f';  f = f(trimN : end - trimN);
% 3. Notch filtering
tw = tukeywin(length(f), 1000/length(f));
f_tw = f.*tw';
f_notch  = filter(F_Notch, f_tw);
% 4. highpass filtering
hp_tw = tukeywin(length(f_notch), 1000/length(f_notch));
f_notch_tw = f_notch.*hp_tw';
f_notch_hp = filter(F_HP, f_notch);

d = {y_md f f_notch f_notch_hp};

figure 
for i = 1:length(d)
%     [sf, xf] = pwelch(d{i}, fs/res, 0.1*fs/res, [], fs);
%     [sf, xf] = pburg( d{i}, 256, Fspan, fs);
    subplot(2,1,1)
    hold on;
    plot(d{i})
    subplot(2,1,2)
    hold on;
    semilogy(abs(fft(d{i})))
end
xlim([-1,60])
subplot(2,1,1)
legend('raw', '+ denoise', '+ notch', '+ hp')

function [ft_dn_notch_hp_rsmp, m_artif1, m_artif2, m_cost, ft_dn, ft_dn_notch, ft_dn_notch_hp]...
    = ft_preproc_sv(dat, highpassB, notchB, Nit, LABEL, svpath, TARApath)
    tic
    % highpassB = 0.1; notchB = 50; Nit = 50;
    [ft_dn_notch_hp_rsmp, m_artif1, m_artif2, m_cost, ft_dn, ft_dn_notch, ft_dn_notch_hp] = ...
        ft_TARA_notch_hp_resample_l(dat, Nit, highpassB, notchB, TARApath);
    save(strcat(svpath, 'ft_dn_notch_hp_rsmp_', LABEL), 'ft_dn_notch_hp_rsmp');
    save(strcat(svpath, 'm_artif1_',            LABEL), 'm_artif1');
    save(strcat(svpath, 'm_artif2_',            LABEL), 'm_artif2');
    save(strcat(svpath, 'm_cost_',              LABEL), 'm_cost');
    save(strcat(svpath, 'ft_dn_',               LABEL), 'ft_dn');
    save(strcat(svpath, 'ft_dn_notch_',         LABEL), 'ft_dn_notch');
    save(strcat(svpath, 'ft_dn_notch_hp_',      LABEL), 'ft_dn_notch_hp');
    toc
end

function [res_dat, artif1, artif2, costM, dn_dat, dn_notch_dat, dn_notch_hp_dat] = ...
        ft_TARA_notch_hp_resample_l(ft_dat, Nit, highpassB, notchB, TARApath )
%% TARA denosing the trial data of a fieldtrip data structure. 
% 2017-11-8 by Vincent
% ft_dat : data with fieldtrip format
% Nit    : number of iteration, usually 50 to 100;
addpath(TARApath);
res_dat         = ft_dat; 
dn_dat          = ft_dat; 
dn_notch_dat    = ft_dat; 
dn_notch_hp_dat = ft_dat;

N_chans = ft_dat.hdr.nChans;
[~, n_sig]  = size(ft_dat.trial{1});
fs = ft_dat.fsample;
trimN = 3;
n_sig_f  = n_sig - 2*trimN + 1;

% fc : cut-off frequency for TARA highpass filter (cycles/sample) (0 < fc < 0.5)
% 0.03 for 4096Hz is 122.88, and 0.12288 for 1000Hz to be 122.88;
if fs == 1000
    RESMP = 0;
    fc = 0.12288;
elseif fs == 4096
    disp('Data need to be resampled to 1000 Hz...')
    RESMP = 1;
    fc = 0.03;    
    y  = ft_dat.trial{1}(1,:);
    yr = resample(y(trimN:end-trimN), 1000, fs);
    n_resmp = length(yr);
    dn_notch_hp_rsmp = zeros(N_chans, n_resmp);
else
    disp('sampling frequency error!')
end

% configuration for TARA 
d = 2;                      % d : filter order parameter (d = 1, 2, or 3)
theta = 0.05; beta  = 1.4;  % TARA initialize parameter 
ps    = 100;                 % psuedo sigma (psuedo noise)

% container
artif1      = zeros(N_chans, n_sig);
artif2      = zeros(N_chans, n_sig);
costM       = zeros(N_chans, Nit);
dn          = zeros(N_chans, n_sig_f);
dn_notch    = zeros(N_chans, n_sig_f);
dn_notch_hp = zeros(N_chans, n_sig_f);

% notch and high pass filter design
% fs =4096; highpassB = 0.1; notchB = 50; Nit = 50; trimN = 100;
F_Notch = designfilt('bandstopiir','FilterOrder',10, ...
         'HalfPowerFrequency1',notchB-2,'HalfPowerFrequency2',notchB+2, ...
         'SampleRate',fs);
f_hp = highpassB/(fs/2);
dHP  = fdesign.highpass('Fst,Fp,Ast,Ap', f_hp*1, f_hp*1.5, 60, 2);
F_HP = design(dHP,'butter');

%% denoise - notch - high pass - resample
for i_ch = (1: N_chans-1)
    i_ch
    y = ft_dat.trial{1}(i_ch,:);
    % 1. de-median
    y_dm  = y - median(y);
    % 2. TARA denosing
    [x1, x2, f, cost] = tara2_L1(y_dm, d, fc, theta, beta, ps, Nit);
    f = f';  f = f(trimN : end - trimN);
    dn(i_ch,:)      = f;
    artif1(i_ch,:)  = x1;
    artif2(i_ch,:)  = x2;
    costM(i_ch,:)   = cost;
    % 3. Notch filtering
    f_notch  = filter(F_Notch, f);
    dn_notch(i_ch,:) = f_notch;
    % 4. highpass filtering
    f_notch_hp = filter(F_HP, f_notch);
    dn_notch_hp(i_ch,:)  = f_notch_hp;
    % 5. resample
    if RESMP == 1
        f_notch_hp_rsmp          = resample(f_notch_hp, 1000, fs);
        dn_notch_hp_rsmp(i_ch,:) = f_notch_hp_rsmp;
    end
end
disp('Redormating results...')
% for demedian and denoise data
dn_dat.trial = {dn}; dn_dat.fsample = 1000; 
trl = []; trl(1) = 1; trl(2) = length(dn_dat.trial{1}); trl(3) = 0;
dn_dat.time = {(1:size(dn_dat.trial{1},2))/dn_dat.fsample};
cfg = []; cfg.trl = trl; dn_dat = ft_redefinetrial(cfg, dn_dat);
disp('Denoised data obtained...')
% for denoise + notch
dn_notch_dat.trial    = {dn_notch};    dn_notch_dat.fsample    = 1000;     
trl = []; trl(1) = 1; trl(2) = length(dn_notch_dat.trial{1}); trl(3) = 0;
dn_notch_dat.time = {(1:size(dn_notch_dat.trial{1},2))/dn_notch_dat.fsample};
cfg = []; cfg.trl = trl; dn_notch_dat = ft_redefinetrial(cfg, dn_notch_dat);
disp('Denoised + notched data obtained...')
% for denoise + notch + high pass
dn_notch_hp_dat.trial = {dn_notch_hp}; dn_notch_hp_dat.fsample = 1000; 
trl = []; trl(1) = 1; trl(2) = length(dn_notch_hp_dat.trial{1}); trl(3) = 0;
dn_notch_hp_dat.time = {(1:size(dn_notch_hp_dat.trial{1},2))/dn_notch_hp_dat.fsample};
cfg = []; cfg.trl = trl; dn_notch_hp_dat = ft_redefinetrial(cfg, dn_notch_hp_dat);
disp('Denoised + notched + highpass data obtained...')
% for denoise + notch + high pass + resample
if RESMP == 1
    res_dat.fsample = 1000; res_dat.trial = {dn_notch_hp_rsmp} ;
    % refresh ft data structur
    trl = []; trl(1) = 1; trl(2) = length(res_dat.trial{1}); trl(3) = 0;
    res_dat.time = {(1:size(res_dat.trial{1},2))/res_dat.fsample};
    cfg = []; cfg.trl = trl; res_dat = ft_redefinetrial(cfg, res_dat);
else
    res_dat.fsample = 1000;
    res_dat.trial = {dn_notch_hp};
    % refresh ft data structur
    trl = [];
    trl(1) = 1; trl(2) = length(res_dat.trial{1}); trl(3) = 0;
    res_dat.time = {(1:size(res_dat.trial{1},2))/res_dat.fsample};
    cfg = []; cfg.trl = trl;
    res_dat = ft_redefinetrial(cfg, res_dat);
end
end
%% prepare lead field
% load('Cortex-mid_Su_reduced');  cortex_mid = nP;   clear nP;
% eeg_LF_r = load('EEG_leadfield_reduced'); ecog_LF_r = load('EcoG_leadfield_reduced'); 
% load('EEG-leadfield'); eeg_LF = LF; clear LF;
% load('EcoG-leadfield'); ecog_LF = LFc; clear LFc;
% eeg_LF=eeg_LF(:,eeg_LF_r.ind'); ecog_LF=ecog_LF(:,ecog_LF_r.ind');
% deriv_order = [1     2     4     6     9    10    13    15 ...
%     17    18     3     7     8    11    12    16     5    14];
% eeg_LF = eeg_LF(deriv_order,:);
% save('EEG_LF.mat', 'eeg_LF')
% save('ECoG_LF.mat', 'ecog_LF')