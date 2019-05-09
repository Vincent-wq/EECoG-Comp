% EECoG_Comp_Sim_min_test_Monkey.m 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EECoG-Comp minimum test bench with simulated simultaneous EEG and ECoG data.
% xxx
% created by Vinceng and Pedro.
% Last modified in 2019.1.12
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear; close all; 
P= cd; FolderLength = length('OHBM-2019');
addpath(genpath([P(1:end-FolderLength),'\tools\QUIC\']));

addpath(genpath([P(1:end-FolderLength),'\dataWD\']));
svPath = [P(1:end-FolderLength),'dataWD\'];
clear P FolderLength;
%% 1.load data
overWrite = 1;
tNow = char(datetime('now','Format','_yyyy_MM_dd'))
% tNow = '_2019_01_08';
resL = '1600'; resH = '10007';
load([svPath 'INFO-' resH '-' resL]); simRes = '_t1_15_128';
% load([svPath 'Sim_' 'INFO_' resH '_' resL '_2019_01_08']);
load([svPath 'Geo-' resH '-' resL ],'ScECoG2sECoG', 'L','ScEEG2sECoG', 'cortexL', 'cortexH', 'indsL');
%% 2. Generate sim data
nSample = 10800; GENS = [150,200,300,350]; nGENS = length(GENS);
 Xeeg=cell(nGENS,1); Yeeg=cell(nGENS,1); Teeg=cell(nGENS,1); AUCeeg=cell(nGENS,1);
Xecog=cell(nGENS,1);Yecog=cell(nGENS,1);Tecog=cell(nGENS,1);AUCecog=cell(nGENS,1);
for xxx = 1:nGENS
    nSour = GENS(xxx);
    PREFIX = ['MiniBench_Monkey_' num2str(nSour) '_' num2str(nSample)];
    info.nSample = nSample; info.nsEEG = nSour;
    para.nTimePoint = info.nSample; para.density = 0.05 ;
    [LFK{1},~,indEEG]   = redLFs(LF{1}{1}, 0,nSour);
    [LFK{2}, ~, indECoG]= redLFs(LF{2}{1}, 0,nSour);
    elePosEEG = LF{1}{3}'; elePosECoG = LF{2}{3}';
    para.eegLF = LFK{1}; para.ecogLF = LFK{2};
    para.eegSNR = 20;  para.ecogSNR= 30; para.senNoise= 1; % 10 20
    para.isComplex = 0; 'CSCG'; para.Show =0;
    para.mode = 'rand'; % 'PMgen' 'CSCG'
    simName = {'Sim EEG'; 'Sim ECoG'};
    RRIM = {'MNE', 'LCMV', 'eLoreta'}'; BLIM = {'eNet-SSBL'}';
    IM = [RRIM; BLIM];
    nRRIM = length(RRIM); nBLIM = length(BLIM); nIM = length(IM);
    % [EEG, ECoG, sData, PM, COV, nz] = EECoG_Comp_Sim_Generator(para,[svPath, 'Sim_data' tNow]);
    [EEG, ECoG, sData, PM, COV, eCOV, nz] = EECoG_Comp_Sim_Generator(para,[]);
    conf = struct('colormap',{'parula'},...
        'ifFig',{1},...
        'linestyle',{'none'},...
        'LineWidth',{4},...
        'alpha',{0.6},...
        'TITLE',{'test Image'},...
        'log',{1},...
        'material',{[0.6, 0.5, 0.3, 1, 0.6]});
    conf.TITLE='Source Space Simulation for EEG and ECoG';
    patchVis(cortexL, cortexL.Vertices(indsL(indEEG),:), diag(PM), conf);
    hold on
    scatter3(elePosEEG(:,1),elePosEEG(:,2),elePosEEG(:,3),50, 'r')
    scatter3(elePosECoG(:,1),elePosECoG(:,2),elePosECoG(:,3),30, 'xm')
    scatter3(elePosECoG(:,1),elePosECoG(:,2),elePosECoG(:,3),30, 'dm')
    testPlot2mat(eCOV, PM,'emp COV', 'true PM', [-0.01 0.01 ], [-0.1 0.1 ], 1,[]);
    
    %  3.2 prepare data
    [ncEEG, nSour]=size(LFK{1}); [ncECoG, ~]=size(LFK{2});
    EEG   = refAve(ncEEG)*EEG;   EEG   = EEG(1:end-1,:);     % average reference and remove one channel
    ECoG  = refAve(ncECoG)*ECoG; ECoG  = ECoG(1:end-1,:);    % average reference and remove one channel
    LFK{1}= refAve(ncEEG)*LFK{1};  LFK{1}=LFK{1}(1:end-1,:); % average reference EEG leadfield
    LFK{2}= refAve(ncECoG)*LFK{2}; LFK{2}=LFK{2}(1:end-1,:); % average reference ECoG leadfield
    dataRes = {EEG;ECoG};  nDataset = length(dataRes);
    info.ncEEG = info.ncEEG-1; info.ncECoG = info.ncECoG-1;
    %% 3. Analysis
    %  3.1 ECoG Laplacian
    [lambda, LAMBDA, GCV] = LambdaGCV(ScECoG2sECoG(1:end-1,:), L, dataRes{2}, 1000, -6, -2, [PREFIX 'ECoG_Laplacian' tNow '_lambda' ]);
    [simH] = calSourLap(ScECoG2sECoG(1:end-1,:), L, lambda);
    simECoGA= (ScEEG2sECoG(indECoG,:)*simH)*ECoG; % simECoGA = simECoGA(1:nSour,:);
    dataECoGA = {simECoGA};  nECoGA = length(dataECoGA);
    % 3.2 solving inverse problem
    DSTF = cell(nDataset,nIM);
    bestLambda = cell(nDataset,nRRIM);
    LAMBDAs = cell(nDataset,nRRIM); GCVs = cell(nDataset,nRRIM);
    lambdaRange = [-3.99, 4]; nLambda = 1000;
    % for MNE, LCMV, and eLoreta
    for iDataset = 1:nDataset
        disp(['Solving inverse problem for ', simName{iDataset}, ' ...'])
        for iIM = 1:nRRIM
            tic
            disp(IM(iIM))
            [GCVs{iDataset, iIM}, LAMBDAs{iDataset, iIM}]= gcvFilter(LFK{iDataset},...
                dataRes{iDataset},RRIM{iIM},lambdaRange, nLambda);
            bestLambda{iDataset, iIM} = LAMBDAs{iDataset,iIM}(GCVs{iDataset, iIM}==min(GCVs{iDataset, iIM}));
            DSTF{iDataset, iIM} = mk_imFilters(dataRes{iDataset},LFK{iDataset}, bestLambda{iDataset, iIM}(1), RRIM{iIM});
            toc
        end
    end
    checkplotIM(LAMBDAs, GCVs, simName, RRIM,1);
    % load([svPath 'Sim_EECoGLeft788_DSTF_2019_01_08']);
%     LIKELIhood = cell(nDataset,nBLIM);
%     PMjj = cell(1,nBLIM);
%     for iDataset = 1:nDataset
%         disp(['Solving inverse problem for ', simName{iDataset}, ' ...'])
%         for iIM = 1:nBLIM
%             tic
%             disp(BLIM(iIM))
%             [tDSTF,tThetaJJ,LIKELIhood{iDataset,iIM}] = mkfilt_iterIM(dataRes{iDataset}, LFK{iDataset}, BLIM{iIM});
%             DSTF{iDataset,nRRIM+iIM} = tDSTF;
%             if iIM == 2; PMjj{1,iIM} = real(tThetaJJ); end
%             toc
%         end
%     end
    % load([svPath 'Sim_bcvareta_directPM' tNow], 'PMjj');
    % %  3.3 alpha band pass
    SimSourEECoG =cell(nDataset, nIM);
    SimSourEECoGAlpha=cell(nDataset, nIM);
    SimLapECoGAlpha = cell(nECoGA,1);
    alpha = [8 12];
    orderFilter = 500;
    disp('Filtering alpha band for Sim EEG, ECoG data...')
    tic
    for iDataset = 1:nDataset
        for iIM = 1:nIM-1
            tmp= DSTF{iDataset, iIM}*dataRes{iDataset};
            SimSourEECoG{iDataset,iIM} = tmp;
            tmp = bpFilter(tmp, info.fs, alpha(1), alpha(2), orderFilter);
            SimSourEECoGAlpha{iDataset, iIM} = tmp;
        end
    end
    toc
    disp('Filtering alpha band for Sim ECoG Laplacian data...')
    tic
    for iECoGA = 1:nECoGA
        SimLapECoGAlpha{iECoGA} = bpFilter(dataECoGA{iECoGA}, info.fs,...
            alpha(1), alpha(2), orderFilter);
    end
    toc
    sDataAlpha = bpFilter(sData, info.fs, alpha(1), alpha(2), orderFilter);
    % %  3.4 PM estimate
    maxIterQUIC = 100;  maxIterSGGM = 60;
    SimEECoGAlphaPM = cell(nDataset, nIM);  SimEECoGAlphaPMn = cell(nDataset, nIM);
    SimLapECoGAlphaPM = cell(nECoGA,1);   SimLapECoGAlphaPMn = cell(nECoGA,1);
    % using Deriel method sggm_lqa
%     SimEECoGAlphaPM_d = cell(nDataset, nIM);  SimEECoGAlphaPM_dn = cell(nDataset, nIM);
%     SimLapECoGAlphaPM_d = cell(nECoGA,1);   SimLapECoGAlphaPM_dn = cell(nECoGA,1);
    
    for iData = 1:nDataset
        for iIM = 1:nIM-1
            disp(['Alpha PM estimation for data ' num2str(iData)...
                num2str(iIM) ' of ' num2str(nDataset*nIM) ' ...'])
            TMP = PM_QUIC(SimSourEECoGAlpha{iData,iIM}, maxIterQUIC);
            SimEECoGAlphaPM{iData,iIM} = TMP;
            SimEECoGAlphaPMn{iData,iIM} = normMat3(TMP);
%             SimEECoGAlphaPM_d{iData,iIM} = real(PMcalVincent(cov(SimSourEECoGAlpha{iData,iIM}'),...
%                 nSample, maxIterSGGM));
%             SimEECoGAlphaPM_dn{iData,iIM} = normMat3(SimEECoGAlphaPM_d{iData,iIM});
        end
    end
    for iECoGA = 1:nECoGA
        disp(['Alpha PM estimation for ECoG Laplacian ' num2str(iECoGA) ' ...'])
        TMP = PM_QUIC(SimLapECoGAlpha{iECoGA,1}, maxIterQUIC);
        SimLapECoGAlphaPM{iECoGA,1}  = TMP;
        SimLapECoGAlphaPMn{iECoGA,1} = normMat3(TMP);
%         SimLapECoGAlphaPM_d{iECoGA,1} = real(PMcalVincent(cov(SimLapECoGAlpha{iECoGA,1}'),...
%             nSample, maxIterSGGM));
%         SimLapECoGAlphaPM_dn{iECoGA,1} = normMat3(SimLapECoGAlphaPM_d{iECoGA,1});
    end
    % source all frequncy pm
    nDataset = 2; nIM = 4;
    for iData = 1:nDataset
        for iIM = 1:nIM-1
            disp(['PM estimation for data ' num2str(iData)...
                num2str(iIM) ' of ' num2str(nDataset*nIM) ' ...'])
            TMP = PM_QUIC(SimSourEECoG{iData,iIM}, maxIterQUIC);
            SimEECoGPM{iData,iIM} = TMP;
            SimEECoGPMn{iData,iIM} = normMat3(TMP);
        end
    end
    for iECoGA = 1
        disp(['Alpha PM estimation for ECoG Laplacian ' num2str(iECoGA) ' ...'])
        TMP = PM_QUIC(dataECoGA{iECoGA,1}, maxIterQUIC);
        SimLapECoGPM{iECoGA,1}  = TMP;
        SimLapECoGPMn{iECoGA,1} = normMat3(TMP);
    end

%     PM_Alpha_sData_QUIC = PM_QUIC(sDataAlpha, maxIterQUIC);
%     PM_Alpha_sData_QUICn = normMat3(PM_Alpha_sData_QUIC);
%     PM_sData_QUIC = PM_QUIC(sData, maxIterQUIC); 
%     PM_sData_QUICn = normMat3(PM_sData_QUIC);
    senPM = cell(2,1); senPMn = cell(2,1);
    for i = 1:2
        TMP = PM_QUIC(dataRes{i,1}, maxIterQUIC);
        senPM{i,1}  = TMP;
        senPMn{i,1} = normMat3(TMP);
    end
    tic
    nBoot = [];
    [Xeeg{xxx},Yeeg{xxx},Teeg{xxx},AUCeeg{xxx}] = calPerfcurve([SimEECoGAlphaPMn(1,1:3), SimLapECoGAlphaPMn],PM,nSample,nBoot);
    toc
    plotPerfcurve(Xeeg{xxx},Yeeg{xxx},AUCeeg{xxx},[IM(1:3);'ECoG LapLacian'], ['Monkey EEG:' ' $n_{c}$=' num2str(ncEEG)...
        ', $n_{s}$=' num2str(nSour) ', $n_{t}$=' num2str(nSample)],...
        ['fig\EECoG_Comp_Sim\\MonkeySimEEG_' para.mode '_' num2str(GENS(xxx)) '_' num2str(nBoot)]);
    
    tic
    [Xecog{xxx},Yecog{xxx},Tecog{xxx},AUCecog{xxx}] = calPerfcurve([SimEECoGAlphaPMn(2,1:3), SimLapECoGAlphaPMn],PM,nSample,nBoot);
    toc
    plotPerfcurve(Xecog{xxx},Yecog{xxx},AUCecog{xxx}, [IM(1:3); {'ECoG LapLacian'}], ['Monkey ECoG:' ' $n_{c}$=' num2str(ncECoG)...
        ', $n_{s}$=' num2str(nSour) ', $n_{t}$=' num2str(nSample)],...
        ['fig\EECoG_Comp_Sim\\MonkeySimECoG_' para.mode '_' num2str(GENS(xxx)) '_' num2str(nBoot)]);
    
    if overWrite; save(['res\MonkeySim_noisecorrected_t' num2str(GENS(xxx)) '_' para.mode simRes tNow]); end
end
if overWrite; save(['res\MonkeySim_resROC_' para.mode tNow],'Xeeg','Yeeg','Teeg','AUCeeg',...
        'Xecog','Yecog','Tecog','AUCecog','GENS', 'nGENS'); end