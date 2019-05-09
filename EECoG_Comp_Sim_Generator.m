% EECoG_Comp_Sim.m 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate simulation data for simultaneous EEG and ECoG comparison. 
% platform.
% created by Vinceng and Pedro.
% Last modified in 2019.1.7
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [EEG, ECoG, sData, PM, COV, eCOV, nz] = EECoG_Comp_Sim_Generator(para, svPath)

n = para.nTimePoint; density = para.density/2;
if isfield(para,'ecogLF') && isfield(para,'ecogSNR'); outECoG =1; else outECoG =0; end
if isfield(para,'senNoise') && para.senNoise ==1; senNoise =1; else senNoise =0; end
eegLF = para.eegLF; 
eegSNR  = db2mag(para.eegSNR); 
if outECoG
    ecogLF = para.ecogLF; 
    ecogSNR = db2mag(para.ecogSNR);
    [nChanECoG, ~]   = size(ecogLF); 
else
    ECoG=[];
end

[nChanEEG, nSour]= size(eegLF); 

if isfield(para,'isComplex') && para.isComplex==1; isComplex = 1; else, isComplex =0; end
if isfield(para,'Show') && para.Show==1; SHOW = 1; else, SHOW =0; end
if isfield(para,'PMclim');   isPMclim = 1; else, isPMclim =0;  end
if isfield(para,'COVclim'); isCOVclim = 1; else, isCOVclim =0; end
% Generate source PM, COV, and data
switch para.mode
    case 'rand'
        [PM, COV, sData, nz] = PMgen(nSour, density, [], n, isComplex, 0.02, 0);
        eCOV = double(cov(sData'));
        EEG  =   eegLF*sData; 
    case 'block'
        if isfield(para,'nBlocks'), nb = para.nBlocks; else, nb = 5; end
        [eCOV, sData, PM] = CSCG(n,nSour,nb,2);
        if isComplex,1; else; PM =real(PM); sData = real(sData)'; eCOV = real(eCOV); end
        COV=speye(nSour)/PM; 
        disp('Sparsity for Cov and PM are: ')
%         figure,plot([std(EEG')',std(EEG1')'])
%         ylim([0, 24])
        nz = [nnz(COV)/(nSour*nSour) nnz(PM)/(nSour*nSour)]
    case 'patch'
        [eCOV, sData, EEG, PM] = preSimPatch({eegLF}, n, para.patchInds, para.indS);
        if isComplex,1; else; PM =real(PM); sData = real(sData)'; 
            EEG=real(EEG);eCOV = real(eCOV); end
        COV=speye(nSour)/PM; 
        disp('Sparsity for patch sim Cov and PM are: ')
        nz = [nnz(COV)/(nSour*nSour) nnz(PM)/(nSour*nSour)]
end
EEG  =   eegLF*sData; 
if senNoise
    AMP = mean(std(EEG'));
    EEG = EEG+AMP*(1/eegSNR)*randn(nChanEEG,n);
end
if outECoG
    ECoG = ecogLF*sData;   
    if senNoise 
        AMP1 = mean(std(ECoG'));
        ECoG = ECoG+AMP1*(1/ecogSNR)*randn(nChanECoG,n);
    end
    %ECoG = ecogLF*sData+(1/ecogSNR)*randn(nChanECoG,n); 
else
    ECoG =[];
end
% plot data
if SHOW
    if isComplex;res1 = abs(COV); res2 = abs(PM);else;res1 = COV; res2 = PM;end
    f1 = figure;
    subplot(221);imagesc(res1);colorbar;if isCOVclim, caxis(para.COVclim);end
    title('COV'); ylabel('Complex Sparse'); 
    subplot(222);imagesc(res2);colorbar;if isPMclim, caxis(para.PMclim); end
    title('PM')
    subplot(223);[res, ~] = preTriHist(res1,1); histogram(res, 100);
    ylabel('Distribution'); title(['nnz rate= ' num2str(nz(1))]);
    subplot(224);[res, ~] = preTriHist(res2,1);histogram(res, 100);
    title(['nnz rate= ' num2str(nz(2))]);
end
% save data
if ischar(svPath)
    set(f1,'Unit','normalized','Position',[0.05,0.05, 6, 1]);
    print(['Sim_COV_' num2str(round(nz(1),4)) '_PM_' num2str(round(nz(2),4)) '.tiff'],'-dtiff','-r300');
    save(svPath, 'PM', 'COV', 'sData', 'nz', 'EEG', 'ECoG');
end
end