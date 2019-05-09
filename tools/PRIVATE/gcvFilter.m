function [GCV, LAMBDAs] = gcvFilter(K,data,FilterType, lambdaRange, nLambda)
%***************************************************************************
% finds optimal regularization parameter lambda by GCV from Wahba (1979)
%      opt_lambda = min G(l)= || P(l) Y ||^2 ./ Trace(n P(l))^2
% Input:            K: leadfield;
%                data: recordings; 
% Ouput:  LAMBDAsdata: lambda vector;
%                 GCV: GCV vector;
%***************************************************************************
[n,~] = size(data);    
GCV = zeros(nLambda,1);
LAMBDAs = logspace(lambdaRange(1),lambdaRange(2),nLambda);
for i=1:nLambda
    %***********************************************************
    %*     GCV(l) = || P(l) Y ||^2 ./ Trace(n P(l))^2
    %***********************************************************
    switch FilterType
        case 'eLoreta'
            tmp = LAMBDAs(i)^2;
            [mkf,~] = eLoretaFilter(K, n*tmp, 500, .0001);
        case 'LCMV'
            [mkf]=mkfilt_lcmv(K,double(corr(data')), LAMBDAs(i)); 
            mkf = mkf';
        case 'MNE'
            [mkf]=mkfilt_minnorm(K,LAMBDAs(i)); 
            mkf = mkf';
    end
    P = (eye(n) - K*mkf);
    res = P*data;
    GCV(i) = n*norm(res,2)^2/(trace(P))^2; % modified according to Wahba definition.
end
end