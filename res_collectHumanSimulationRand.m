% summarize results for the human random simulations.
clc; clear; close all;
P= cd; FolderLength = length('OHBM-2019');
svPath = [P(1:end-FolderLength),'dataWD\'];

load([P '\res\HumanSim_resROC_rand_2019_01_23'])
invMethods = {'MNE', 'LCMV', 'ELORETA', 'eNet-SSBL'};
AUC =AUCeeg; X =Xeeg; Y =Yeeg;
[nSen, nSour] = size(AUC); nM = length(invMethods);
%% corrected AUCx
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
colours = permute(get(gca, 'colororder'), [1 3 2]);
for i =1:4
    subplot(2,2,i)
    imagesc(squeeze(datMat(i,:,:)),'AlphaData' ,0.7)
    xticklabels({'150';'200';'300';'350'});
    yticklabels({'1020';'32';'64';'128'})
    caxis([0,1])
    switch i
        case 1
            hx=xlabel('MNE','interpreter','latex','FontSize',14);
            v=get(hx,'Position'); v(2)=v(2)+0.1; set(hx,'Position',v);
            hand_t=title('Human EEG ESI Methods Performance for Random Source ($spAUC=0.1$)'...
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
%print('fig\HumanSim_ESIM_random1.tiff', '-dtiff','-r300');



