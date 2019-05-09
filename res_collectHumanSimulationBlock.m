% summarize results for the human block simulations.
clc; clear; close all;
P= cd; FolderLength = length('OHBM-2019');
svPath = [P(1:end-FolderLength),'dataWD\'];

load([P '\res\HumanSenSourTestXYAUC20190120'])
invMethods = {'MNE', 'LCMV', 'ELORETA', 'eNet-SSBL'};
[nSen, nSour] = size(AUC); nM = length(invMethods);
CHANS = [1020,32,64,128]; GENS = [150,200,300,350];

% pos1 = zeros(nSen, nSour); low1 = zeros(nSen, nSour); high1 = zeros(nSen, nSour);
% pos2 = zeros(nSen, nSour); low2 = zeros(nSen, nSour); high2 = zeros(nSen, nSour);
% pos3 = zeros(nSen, nSour); low3 = zeros(nSen, nSour); high3 = zeros(nSen, nSour);
% pos4 = zeros(nSen, nSour); low4 = zeros(nSen, nSour); high4 = zeros(nSen, nSour);
% for iSen = 1:nSen
%     for iSour = 1:nSour
%         ele = AUC{iSen,iSour}; ele1 = ele{1}; ele2 = ele{2}; ele3 = ele{3}; ele4 = ele{4};
%         pos1(iSen,iSour) = ele1(1); low1(iSen,iSour) = pos1(iSen,iSour)-ele1(2); high1(iSen,iSour) = ele1(3)-pos1(iSen,iSour);
%         pos2(iSen,iSour) = ele2(1); low2(iSen,iSour) = pos2(iSen,iSour)-ele2(2); high2(iSen,iSour) = ele2(3)-pos2(iSen,iSour);
%         pos3(iSen,iSour) = ele3(1); low3(iSen,iSour) = pos3(iSen,iSour)-ele3(2); high3(iSen,iSour) = ele3(3)-pos3(iSen,iSour);
%         pos4(iSen,iSour) = ele4(1); low4(iSen,iSour) = pos4(iSen,iSour)-ele4(2); high4(iSen,iSour) = ele4(3)-pos4(iSen,iSour);
%     end
% end
%% new plot across sensors and sources
% [posBase1,posBase2] = meshgrid(0:0.2:0.6, 0:0.2:0.6);
% f=figure('Units','normalized','Position',[0.05,0.05, 0.6, 0.7]);
% colours = permute(get(gca, 'colororder'), [1 3 2]);
% x1=(posBase1+pos1); y1 = (posBase2+pos1); x2=(posBase1+pos2); y2 = (posBase2+pos2);
% x3=(posBase1+pos3); y3 = (posBase2+pos3); x4=(posBase1+pos4); y4 = (posBase2+pos4);
% yStep = posBase2(:,1)+0.5;
% hold on;
% grid on;
% er1=errorbar(x1(:),y1(:),low1(:),high1(:),low1(:),high1(:),'Color',...
%     squeeze(colours(1,:,:))','Marker','o','LineStyle','none');
% er2=errorbar(x2(:),y2(:),low2(:),high2(:),low2(:),high2(:),'Color',...
%     squeeze(colours(2,:,:))','Marker','o','LineStyle','none');
% er3=errorbar(x3(:),y3(:),low3(:),high3(:),low3(:),high3(:),'Color',...
%     squeeze(colours(3,:,:))','Marker','o','LineStyle','none');
% er4=errorbar(x4(:),y4(:),low4(:),high4(:),low4(:),high4(:),'Color',...
%     squeeze(colours(4,:,:))','Marker','o','LineStyle','none');
% legend([er1,er2,er3,er4],invMethods)
% xlim([0.5, 1.4])
% ylim([0.5, 1.3])
% xticks(0.5:0.1:1.3) % [yStep;1.3]+0.1
% xticklabels({'','150','','200','','300','','350'})
% yticks(0.5:0.1:1.3) % yStep+0.1
% yticklabels({'','1020','','32','','64','','128'}) % {num2str(CHANS')}
% ylabel({'$n_{sensor}$, $AUC=0.6\pm0.1$'},'Interpreter','latex')
% xlabel({'$n_{source}$, $AUC=0.6\pm0.1$'},'Interpreter','latex')
% title('AUC for Block Source Settings','Interpreter','latex')
% set(f,'Unit','normalized','Position',[0.05,0.05, 0.6, 0.8]);
% print('fig\HumanSim_ESIM_block_Matrix.tiff', '-dtiff','-r300');

%% corrected AUCx
colours = permute(get(gca, 'colororder'), [1 3 2]);
% plot partial AUC for sensor*source*invMethods
P =0.1;
[XeegP,YeegP,AUCP]    = calPartialROC(X,Y, P);
dataMat1 =[cell2mat(AUCP{1,1})',cell2mat(AUCP{1,2})',cell2mat(AUCP{1,3})',cell2mat(AUCP{1,4})'];
dataMat2 =[cell2mat(AUCP{2,1})',cell2mat(AUCP{2,2})',cell2mat(AUCP{2,3})',cell2mat(AUCP{2,4})'];
dataMat3 =[cell2mat(AUCP{3,1})',cell2mat(AUCP{3,2})',cell2mat(AUCP{3,3})',cell2mat(AUCP{3,4})'];
dataMat4 =[cell2mat(AUCP{4,1})',cell2mat(AUCP{4,2})',cell2mat(AUCP{4,3})',cell2mat(AUCP{4,4})'];
datMat = zeros(nM, nSen, nSour);
for i = 1:4
    datMat(i,:,:) = [dataMat1(i,:);dataMat2(i,:);dataMat3(i,:);dataMat4(i,:)];
end
figure
for i =1:4
    subplot(2,2,i)
    imagesc(squeeze(datMat(i,:,:)),'AlphaData' ,1)
    xticklabels({'150';'200';'300';'350'});
    yticklabels({'1020';'32';'64';'128'})
    caxis([0.8,1])
    switch i
        case 1
            hx=xlabel('MNE','interpreter','latex','FontSize',14);
            v=get(hx,'Position'); v(2)=v(2)+0.1; set(hx,'Position',v);
            hand_t=title('Human EEG ESI Methods Performance for Block Source ($spAUC=0.1$)'...
                ,'FontSize',16,'interpreter','latex');
            v=get(hand_t,'Position');v(1)=v(1)+2.6; v(2)=v(2)-0.3;
            set(hand_t,'Position',v);
        case 2
             hx=xlabel('LCMV','interpreter','latex','FontSize',14);
             v=get(hx,'Position'); v(2)=v(2)+0.1; set(hx,'Position',v);
        case 3
            hx=xlabel('eLORETA','interpreter','latex','FontSize',14);
            v=get(hx,'Position'); v(2)=v(2)+0.1; set(hx,'Position',v);
        case 4
            hx=xlabel('eNet-SSBL','interpreter','latex','FontSize',14);
            v=get(hx,'Position'); v(2)=v(2)+0.1; set(hx,'Position',v);
    end
end
set(gcf,'Unit','normalized','Position',[0.05,0.05, 0.68, 0.9]);
hp4 = get(subplot(2,2,4),'Position');
colorbar('Position', [hp4(1)+0.37 hp4(2)  0.05  0.82],...
    'Ticks',[0,0.2, 0.4, 0.5,0.6,0.8,1],...
    'TickLabels',{'0';'0.2';'0.4';'0.5';'0.6';'0.8';'1.0'},...
    'AxisLocation','in');
%colormap(bipolar)
colormap(hot)
%print('fig\HumanSim_ESIM-block1.tiff', '-dtiff','-r300');

%% collect variables
% CHANS = [19,32,64,128]; GENS = [150,200,300,350];
% nCHANS = length(CHANS); nGENS = length(GENS);
% Xeeg=cell(nCHANS, nGENS); Yeeg=cell(nCHANS, nGENS);
% Teeg=cell(nCHANS, nGENS); AUCeeg=cell(nCHANS, nGENS);
%
% load([ P '\res\HumanSim_block_t0_19_200_2019_01_20.mat'], 'Xeeg', 'Yeeg', 'AUCeeg');
% AUC = AUCeeg; X = Xeeg; Y = Yeeg;
% load([ P '\res\HumanSim_block_t0_19_300_2019_01_20.mat'], 'Xeeg', 'Yeeg', 'AUCeeg');
% X{1,3}= Xeeg; Y{1,3}= Yeeg; AUC{1,3}= AUCeeg;
% load([ P '\res\HumanSim_block_t0_19_350_2019_01_20.mat'], 'Xeeg', 'Yeeg', 'AUCeeg');
% X{1,4}= Xeeg{1,3}; Y{1,4}= Yeeg{1,3}; AUC{1,4}= AUCeeg{1,3};
%
% load([ P '\res\HumanSim_block_t0_32_200_2019_01_20.mat'], 'Xeeg', 'Yeeg', 'AUCeeg');
% X(2,1:2)= Xeeg(2,1:2); Y(2,1:2)= Yeeg(2,1:2); AUC(2,1:2)= AUCeeg(2,1:2);
% load([ P '\res\HumanSim_block_t0_32_300_2019_01_20.mat'], 'Xeeg', 'Yeeg', 'AUCeeg');
% X{2,3}= Xeeg; Y{2,3}= Yeeg; AUC{2,3}= AUCeeg;
% load([ P '\res\HumanSim_block_t0_32_350_2019_01_20.mat'], 'Xeeg', 'Yeeg', 'AUCeeg');
% X{2,4}= Xeeg{2,3}; Y{2,4}= Yeeg{2,3}; AUC{2,4}= AUCeeg{2,3};
%
% load([ P '\res\HumanSim_block_t0_64_350_2019_01_20.mat'], 'Xeeg', 'Yeeg', 'AUCeeg');
% X(3,1:2)= Xeeg(3,1:2); Y(3,1:2)= Yeeg(3,1:2); AUC(3,1:2)= AUCeeg(3,1:2);
% X{3,4}= Xeeg{3,3}; Y{3,4}= Yeeg{3,3}; AUC{3,4}= AUCeeg{3,3};
% load([ P '\res\HumanSim_block_t0_64_300_2019_01_20.mat'], 'Xeeg', 'Yeeg', 'AUCeeg');
% X{3,3}= Xeeg; Y{3,3}= Yeeg; AUC{3,3}= AUCeeg;
%
% load([ P '\res\HumanSim_block_t0_128_350_2019_01_20.mat'], 'Xeeg', 'Yeeg', 'AUCeeg');
% X(4,1:2)= Xeeg(4,1:2); Y(4,1:2)= Yeeg(4,1:2); AUC(4,1:2)= AUCeeg(4,1:2);
% X{4,4}= Xeeg{4,3}; Y{4,4}= Yeeg{4,3}; AUC{4,4}= AUCeeg{4,3};
% load([ P '\res\HumanSim_block_t0_128_300_2019_01_20.mat'], 'Xeeg', 'Yeeg', 'AUCeeg');
% X{4,3}= Xeeg; Y{4,3}= Yeeg; AUC{4,3}= AUCeeg;
% X(5,:) =[]; X(:,5) =[]; Y(5,:) =[]; Y(:,5) =[]; AUC(5,:) =[]; AUC(:,5) =[];
% save('senSourTestXYAUC20190120', 'X', 'Y', 'AUC','CHANS','GENS')


