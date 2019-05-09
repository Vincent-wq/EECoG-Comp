function [T, lambdaOptimal] = eLoretaFilter(K,lambda,maxIter,tolerance)
% makes spatial filter according to eLoreta
% input:   K : nChan*nSour leadfield for nChan channels, nSour sources
%      lambda: optional regularization parameter;
%     maxIter: Maximum number of iteration;
%   tolerance: the max error to jump out itertation.
% output:  W : nChan*nSour Filter, nChan channels, nSour sources.
% X Eduardo Gonzalez-Moreira, April 2017
[nchan,nSour] = size(K);
W = ones(1,nSour);
H = eye(nchan)-ones(nchan)/nchan;
No_Conv = 1;
%data = [];
ic = 1;
while No_Conv
    W_old = W;
    WK = diag(W)\K';
    lambdaOptimal = lambda*trace((K*WK))/nchan;
    M = pinv((K*WK)+lambdaOptimal*H);
    W = sqrt(mean(K'*M.*K',2))';
    %W_conv = sqrt(sum(abs(W-W_old).^2));
    W_old_norm = W_old./max(abs(W_old(:)));
    W_norm     = W./max(abs(W(:))); 
    W_conv     = sqrt(sum(abs(W_norm-W_old_norm).^2));  %Vincent
    %data(end+1)=W_conv;
    %fprintf('.');
    ic=ic+1;
    if ic == maxIter          %500
        KM = K'*M;
        T = diag(W)\KM; 
        No_Conv = 0;
    end
    if W_conv < tolerance   % .000001
        KM = K'*M;
        T = diag(W)\KM;        
        No_Conv = 0;
    end
end
end