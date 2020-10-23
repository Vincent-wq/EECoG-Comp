function PM = PM_QUIC(data, maxIter, isStand)
[p n]=size(data);
if isStand
    eCov  = double(corr(data'));  % corr for statistical properties
else
    eCov  = double(cov(data'));  % corr for statistical properties
end

lambda    = sqrt(log(p)/n);                 % Regularization parameter                                                    % Regularization mask nondiagonal
lambdaMat = (lambda*(ones(p)-eye(p))).^2;   % Regularization mask squared
lambdaMat = lambdaMat.*(real(diag(eCov))*real(diag(eCov))');
%maxIter = 80;
[PM, Coveeg, opt, time, iter, dGap] = QUIC('default', eCov, lambdaMat, 1e-3, 1, maxIter);
time
%[PM] = debiasPM(PM, eCov);
end 