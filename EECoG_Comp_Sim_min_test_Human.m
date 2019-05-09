% EECoG_Comp_Sim_min_test_Human.m 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EECoG-Comp minimum test bench with Human forward model, take left EEG
% as EEG and Full EEG as ECoG.
% xxx
% created by Vinceng and Pedro.
% Last modified in 2019.1.13
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
load([svPath 'INFO-' resH '-' resL],'info');
%% HuamModel_15_150,  HuamModel_128_150
% simulation run for different number of sensors and sources.
% ,400,450
CHANS = [19,32,64,128]; GENS = [150,200,300,350];
nCHANS = length(CHANS); nGENS = length(GENS);
Xeeg=cell(nCHANS, nGENS); Yeeg=cell(nCHANS, nGENS); 
Teeg=cell(nCHANS, nGENS); AUCeeg=cell(nCHANS, nGENS);
nSample = 10800; 
nSimSourFull = 1604; 
for xxx=1:nCHANS
    nSimChan = CHANS(xxx)
    for yyy = 1:nGENS
        nSimSour = GENS(yyy)
        load(['HumanHeadModel_' num2str(nSimChan) '_' num2str(nSimSourFull)]);
        simRes = ['_t0_' num2str(nSimChan) '_' num2str(nSimSour)];
        %% 2. Generate sim data
        [LF,~,indEEG]= redLFs(LF, nSimChan, nSimSour);
        para.eegLF  = refAve(nSimChan)*LF;  % average reference for leadfield
        [ncEEG,~]=size(LF);
        info.nSample = nSample; info.nsEEG = nSimSour;
        para.nTimePoint = info.nSample; para.density = 0.05 ;
        para.eegSNR = 10; para.patchInds = patchIndsSel; %para.indS = indS;
        para.isComplex = 0; para.Show =0;
        para.mode = 'rand'; %'block' 'PMgen' 'CSCG' 'patch'
        [EEG, ECoG, sData, PM, COV, eCOV, nz] = EECoG_Comp_Sim_Generator(para,[]);
        conf = struct('colormap',{'parula'},...
            'ifFig',{1},...
            'linestyle',{'none'},...
            'LineWidth',{4},...
            'alpha',{0.6},...
            'TITLE',{'test Image'},...
            'log',{1},...
            'material',{[0.6, 0.5, 0.3, 1, 0.6]});
        conf.TITLE=['Source Space Simulation for ' num2str(nSimChan) '-' num2str(nSimSour) ' Human EEG'];
        patchVis(cortexHuman, cortexHuman.Vertices(indEEG,:), diag(PM), conf);
        hold on
        scatter3(elePos(:,1),elePos(:,2),elePos(:,3),50,'r')
        testPlot2mat(COV, PM,'emp COV', 'true PM', [-0.01 0.01 ], [-0.1 0.1 ], 1,[]);
        % take average reference
        simName = {'Human Sim EEG'}; info.EECoGname = simName;
        RRIM = info.RRIM; BLIM = info.BLIM(1); IM = info.invMethods(1:4);
        nRRIM = length(RRIM); nBLIM = length(BLIM); nIM = length(IM);
        
        EEG   = refAve(ncEEG)*EEG;   EEG   = EEG(1:end-1,:);  % average reference and remove one channel
        dataRes = {EEG}; nDataset = length(dataRes);
        LFK = {LF(1:end-1,:)}; info.ncEEG = info.ncEEG-1; % format leadfield
        %% 3. solving inverse problems
        DSTF = cell(nDataset,nRRIM+nBLIM);
        bestLambda = cell(nDataset,nRRIM);
        LAMBDAs = cell(nDataset,nRRIM); GCVs = cell(nDataset,nRRIM);
        lambdaRange = [-4, 2]; nLambda = 1000;
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
        % % load([svPath 'Sim_EECoGLeft788_DSTF_2019_01_08']);
        LIKELIhood = cell(nDataset,nBLIM);
        PMjj = cell(1,nBLIM);
        for iDataset = 1:nDataset
            disp(['Solving inverse problem for ', simName{iDataset}, ' ...'])
            for iIM = 1:nBLIM
                tic
                disp(BLIM(iIM))
                [tDSTF,tThetaJJ,LIKELIhood{iDataset,iIM}] = mkfilt_iterIM(dataRes{iDataset},...
                    LFK{iDataset}, BLIM{iIM});
                DSTF{iDataset,nRRIM+iIM} = tDSTF;
                PMjj{1,iIM} = real(tThetaJJ);
                toc
            end
        end
        % %  3.3 alpha band pass
        SimSourEECoG =cell(nDataset, nIM);
        SimSourEECoGAlpha=cell(nDataset, nIM);
        alpha = [8 12];
        orderFilter = 500;
        disp('Filtering alpha band for Sim ECoG data...')
        tic
        for iDataset = 1:nDataset
            for iIM = 1:nIM
                tmp= DSTF{iDataset, iIM}*dataRes{iDataset};
                SimSourEECoG{iDataset,iIM} = tmp;
                tmp = bpFilter(tmp, info.fs, alpha(1), alpha(2), orderFilter);
                SimSourEECoGAlpha{iDataset, iIM} = tmp;
            end
        end
        toc
        sDataAlpha = bpFilter(sData, info.fs, alpha(1), alpha(2), orderFilter);
        % %  3.4 PM estimate
        maxIterQUIC = 100; maxIterSGGM = 60;
        SimEECoGAlphaPM = cell(nDataset, nIM);  SimEECoGAlphaPMn = cell(nDataset, nIM);
        % using Deriel method sggm_lqa
        SimEECoGAlphaPM_d = cell(nDataset, nIM);  SimEECoGAlphaPM_dn = cell(nDataset, nIM);
        for iData = 1:nDataset
            for iIM = 1:nIM
                disp(['Alpha PM estimation for data ' num2str(iData) ' with Inverse Solution '...
                    num2str(iIM) ' of ' num2str(nDataset*nIM) ' ...'])
                TMP = PM_QUIC(SimSourEECoGAlpha{iData,iIM}, maxIterQUIC);
                SimEECoGAlphaPM{iData,iIM} = TMP;
                SimEECoGAlphaPMn{iData,iIM} = normMat3(TMP);
                SimEECoGAlphaPM_d{iData,iIM} = real(PMcalVincent(cov(SimSourEECoGAlpha{iData,iIM}'),...
                    nSample, maxIterSGGM));
                SimEECoGAlphaPM_dn{iData,iIM} = normMat3(SimEECoGAlphaPM_d{iData,iIM});
            end
        end
        sPM = PM_QUIC(sDataAlpha, maxIterQUIC); sPMs = normMat3(sPM);
        % PM_Alpha_sData_QUIC = PM_QUIC(sDataAlpha, maxIterQUIC);
        % PM_Alpha_sData_QUICn = normMat3(PM_Alpha_sData_QUIC);
        % PM_Alpha_sData_SGGM = real(PMcalVincent(cov(sDataAlpha'), nSample, maxIterSGGM));
        % PM_Alpha_sData_SGGMn = normMat3(PM_Alpha_sData_SGGM);
        %
        % PM_sData_QUIC = PM_QUIC(sData, maxIterQUIC); PM_sData_QUICn = normMat3(PM_sData_QUIC);
        % PM_sData_SGGM = real(PMcalVincent(cov(sData'), nSample, maxIterSGGM));
        % PM_sData_SGGMn = normMat3(PM_sData_SGGM);
        
        %% final test
        tic
        nBoot = []; %200
        [Xeeg{xxx,yyy},Yeeg{xxx,yyy},Teeg{xxx,yyy},AUCeeg{xxx,yyy}] = calPerfcurve(SimEECoGAlphaPMn,sPMs,nSample,nBoot);
        toc
        plotPerfcurve(Xeeg{xxx,yyy},Yeeg{xxx,yyy},AUCeeg{xxx,yyy},IM, ['Human EEG:' '$n_{c}$=' num2str(nSimChan)...
            ', $n_{s}$=' num2str(nSimSour) ', $n_{t}$=' num2str(nSample)],...
            ['fig\EECoG_Comp_Sim\\HumanSim_' para.mode '_' simRes(5:end) '_' num2str(nBoot)]);
        
        if overWrite; save(['res\HumanSim_' para.mode simRes tNow]); end
    end
end
if overWrite; save(['res\HumanSim_resROC_' para.mode tNow],'Xeeg','Yeeg','Teeg','AUCeeg',...
        'CHANS','GENS','tNow','nCHANS','nGENS','info','nSample','nSimSourFull'); end
