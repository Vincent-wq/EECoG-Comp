function [PMn] = normMat3(PM)
% normalize a covariance or precision matrix by the square of diags.
% by Vincent @ 10th, Dec. 2018 Cuba
[nChan,~,nF]=size(PM);
PMn = zeros(nChan, nChan, nF);
for k = 1:nF
    tmp =squeeze(PM(:,:,k));
    tmpD = diag(tmp);
    tmpDD = diag(sqrt(tmpD));
    PMn(:,:,k) = tmpDD\tmp/tmpDD;
end
end