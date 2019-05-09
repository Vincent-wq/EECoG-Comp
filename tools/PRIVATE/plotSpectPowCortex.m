function [] = plotSpectPowCortex(F, Svv, H, cortex, inds, VIDEO, NAME, conf)
%% check laplacian pow spect distribution for all freqs
fL = length(F);
[nsL,~] = size(H);
if VIDEO
    aviobj = VideoWriter([NAME,'-powSpect.avi']);
    open(aviobj);
end
mkdir(NAME);
for i =1:fL
    COV = squeeze(Svv(:,:,i));
    [covsq] = spd2pd(sqrtm(COV));
    LeftCov = H*covsq;
    RES = zeros(1,nsL);
    for j =1: nsL
        RES(j)=real(LeftCov(j,:)*LeftCov(j,:)');
    end
    [~, indMaxPow] = max(RES)
    figure
    patchVis(cortex, cortex.Vertices(inds,:), RES, conf);
    hold on
    scatter3(cortex.Vertices(inds(indMaxPow),1),cortex.Vertices(inds(indMaxPow),2),cortex.Vertices(inds(indMaxPow),3),75,'b','filled');
    view(-60,20);
    title(['Pow. Distr. of ', num2str(F(i)),' Hz @ source: ', num2str(indMaxPow)]);
    fig=gcf;
    FILEPATH = strcat(NAME,'\', num2str(i),'-LapAT-', num2str(F(i)), 'Hz', '.jpg');
    saveas(fig,FILEPATH);
    if VIDEO
        image = imread(FILEPATH);
        writeVideo(aviobj,image);
    end
    close all;
end
if VIDEO
    close(aviobj);
end
end