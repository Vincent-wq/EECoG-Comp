function [DSTF,ThetaJJ,likelihood] = mkfilt_iterIM(data, K, invMethod)
% Create inverse solution DFTF.
% By Vincent, 2019.1.2
[n,nSample] = size(data);   
Svv = double((corr(data')));
switch invMethod
    case 'eNet-SSBL'
        % TMP = cross_nonovgrouped_enet_ssbl({Svv},{K},nSample,preGroup(n));
        % DSTF = TMP{1};
        [DSTF, likelihood] = screening_ssbl(Svv,K,nSample,preGroup(n));
        ThetaJJ=[];
    case 'BC-VARETA'
       [ThetaJJ,DSTF,likelihood] = bcvareta(Svv,K,nSample); 
end
end