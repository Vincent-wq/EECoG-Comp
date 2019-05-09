function [res, f, SPCTMATRIX]=CalMatixSpecF(data, Fs, FREQ, winSize)
% calculate spectrum using Thompson Windows multi-tapper method, for data
% row by row, and find the specific FREQ power for all the channels.
%       input: data: nChan*nPts;
%              Fs  : sampling frequency;
%              FREQ: target frequncy
%           winSize: the trial size for Robust spectrum estimation
%       output: res: target power for all channels at FREQ
%                 f: frequency points
%        SPCTMATRIX: full power matrix for all frequencies and all
%                    channels, nChan*nPowers
% By Vincent 2018.10.26
%% test
dataRaw = data;
[nChan,nPts]=size(dataRaw);
Nt = winSize*Fs;
Ntrial = floor(nPts/Nt);
dataRaw = dataRaw(:,1:Nt*(Ntrial-1));

% remove last trial
t = (0:Nt*Ntrial-Nt-1)/Fs;
params=struct('tapers',[3 5],'Fs',Fs,'trialave',1,'fpass',[1,40]);
quantile_h = .5;
method = struct('class', 'two-tier','tier', struct('estimator',...
    {'mean', 'quantile'}, 'params', {struct(), struct('h',quantile_h)}));
% spec topo plot
res = zeros(1,nChan);

%figure
for i =1:nChan
    disp(['Channel: ' num2str(i)])
	i
    y = reshape(dataRaw(i,:),[Nt, Ntrial-1]);
    [S,f]=mtspectrumc_Robust(y,params, method);
    if i==1
        f = f;
        SPCTMATRIX=zeros(nChan,length(S));
    end
    SPCTMATRIX(i,:)=S;
    ind = find(abs(f-FREQ)==min(abs(f-FREQ)));
    P = S(ind);
    res(1,i) = P;
    %subplot(5,5,cord1020(i))
    %plot_vector(S,f,'l',[],'b',1.2)
    %ylim([-10,30])
    %title(data.label{i})
    %saveas(gcf,['figs/',meta.chan{i}],'jpg');
    %savefig([meta.chan{i},'.fig'])
end

end