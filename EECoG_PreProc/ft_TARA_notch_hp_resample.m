function [res_dat, artif1, artif2, costM, dn_dat, dn_notch_dat, dn_notch_hp_dat] = ft_TARA_notch_hp_resample(ft_dat, Nit, highpassB, notchB )
%% TARA denosing the trial data of a fieldtrip data structure. 
% 2017-11-8 by Vincent
% ft_dat : data with fieldtrip format
% this function is to clean the data: 
% 1. Remove (sync) Trigger channel;
% 2. de-mean;
% 3. denoising with TARA;
% 4. Notch filter for 50;
% 5. High pass filter 0.3;
% Nit    : number of iteration, usually 50 to 100;
addpath('F:\00TychoMonkey\20170920_TychoMonkey_Vincent\Tools\TARA_software');

% remove synch channel
selchan = ft_channelselection({'all' '-Trigger'}, ft_dat.label);
ft_dat = ft_selectdata(ft_dat, 'channel', selchan);

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
