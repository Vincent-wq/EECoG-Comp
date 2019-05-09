function [GCV, LAMBDAs] = gcv_eLoreta(K,data)
%***************************************************************************
% finds optimal regularization parameter lambda by GCV from Wahba (1979)
%      opt_lambda = min G(l)= || P(l) Y ||^2 ./ Trace(n P(l))^2
% Input:            K: leadfield;
%                data: recordings; 
% Ouput:  LAMBDAsdata: lambda vector;
%                 GCV: GCV vector;
%***************************************************************************
[n,~] = size(data);    
npoints = 1000;          % lambda vector;
GCV = zeros(npoints,1);
LAMBDAs = logspace(-8,3,npoints);
for i=1:npoints
    %***********************************************************
    %*     GCV(l) = || P(l) Y ||^2 ./ Trace(n P(l))^2
    %***********************************************************
    tmp = LAMBDAs(i)^2;
    [mkf,~] = eLoretaFilter(K, n*tmp, 500, .0001);
    P = (eye(n) - K*mkf);
    res = P*data;
    GCV(i) = norm(res,2)^2/(trace(P))^2;
end
end