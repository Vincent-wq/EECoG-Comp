% summarize results for the monkey simulations.
clc; clear; close all;
P= cd; FolderLength = length('OHBM-2019');
svPath = [P(1:end-FolderLength),'dataWD\'];

load([P '\res\MonkeySim_resROC_rand_2019_01_24'])
invMethods = {'MNE', 'LCMV', 'ELORETA', 'ECoG Laplacian'};
nM = length(invMethods);

%% corrected AUCx
colours = permute(get(gca, 'colororder'), [1 3 2]);
% plot partial AUC for sensor*source*invMethods
P =0.1;
[XeegP,YeegP,AUCeegP]    = calPartialROC(Xeeg,Yeeg,P);
[XecogP,YecogP,AUCecogP] = calPartialROC(Xecog,Yecog,P);

dataMat =[cell2mat(AUCeegP{1,1})',cell2mat(AUCeegP{2,1})',...
    cell2mat(AUCeegP{3,1})',cell2mat(AUCeegP{4,1})',...
    cell2mat(AUCecogP{1,1})',cell2mat(AUCecogP{2,1})',...
    cell2mat(AUCecogP{3,1})',cell2mat(AUCecogP{4,1})',];
%% report summary statistics
[min(dataMat(:)) max(dataMat(:))  mean(dataMat(:)) std(dataMat(:))]
disp('Laplacian')
[min(dataMat(4,:)) max(dataMat(4,:))  mean(dataMat(4,:)) std(dataMat(4,:))]
%% group bar plot for different datasets
figure
b = bar(dataMat','FaceAlpha',0.6);
% b.FaceColor =  repmat(squeeze(colours(1:4,:,:)),4,1);
title('Monkey EEG-ECoG ESI Methods Performance for Random Source ($spAUC=0.1$)','interpreter','latex')
xticks(1:1:8)
xticklabels({'EEG 150';'EEG 200';'EEG 300';'EEG 400';'ECoG 150';...
    'ECoG 200';'ECoG 300';'ECoG 400'});
legend(invMethods')
ylim([0.4 0.6]);
xlim([0.5,8.5])
set(gcf,'Unit','normalized','Position',[0.05,0.05, 0.55, 0.8]);
print(['fig\MonkeySimESI_Rand_' num2str(P) '_pAUC.tiff'],'-dtiff','-r300');



