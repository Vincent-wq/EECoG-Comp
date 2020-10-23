function [Cvv] = normCrossSpectrum(Svv)
% input: Svv  : CrossSpectrum: nSensor * nSensor * nFrequncy Points;
% output Cvv  : Normalized Svv by its diagnal.
% by Vincent 2019.4.30
[NSensor,~,Nf] = size(Svv);
Cvv = zeros(NSensor,NSensor,Nf);  % init of the coherence matrix.
for k = 1:Nf
    tmp =squeeze(Svv(:,:,k));
    tmpD = diag(tmp);
    tmpDD = diag(sqrt(tmpD));
    Cvv(:,:,k) = tmpDD\tmp/tmpDD;
end
end