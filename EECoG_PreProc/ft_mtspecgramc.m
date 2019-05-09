function [S,t,f] = ft_mtspecgramc(dat, fs, tapers, fpass, wind, figName, visible)
%% test the multitapper parameters
% fs    : sampling frequency
% tapers: [TW K] where TW is the time-bandwidth product and K is the 
%         number of tapers to be used (less than or equal to 2TW-1).[4, 7]
% fpass:  pass band [0 45]
% wind :  [windSize step]
% figName:save fig path, null for not saving
% visible: visible of figure
% By Vincent 2017.11.24
params.Fs = fs; params.tapers=tapers; 
params.trialave = 0; params.fpass = fpass;
data =  dat.trial{1};
data = data';
[S,t,f] = mtspecgramc(data, wind, params );
tp = size(dat.trial{1},1);
if tp == 19
    S1 = S(:,:,1); % Fp1
    S8 = S(:,:,8); % T3
    S13 = S(:,:,13); %P3
    S17 = S(:,:,17); % O1
elseif tp == 129
    S1 = S(:,:,10); % Fp1
    S8 = S(:,:,69); % T3
    S13 = S(:,:,88); %P3
    S17 = S(:,:,120); % O1
end
h = figure('visible',visible);
subplot(221)
plot_matrix(S1,t,f,'l')
caxis([0,40])
title('Fp1')
subplot(222)
plot_matrix(S8,t,f,'l')
caxis([0,40])
title('T3')
subplot(223)
plot_matrix(S13,t,f,'l')
caxis([0,40])
title('P3')
subplot(224)
plot_matrix(S17,t,f,'l')
caxis([0,40])
title('O1')
if figName ~= ''
    saveas(h, figName)
end