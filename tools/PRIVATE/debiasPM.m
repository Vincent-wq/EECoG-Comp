function [ThetaJJ] = debiasPM(ThetaJJ, Svv)
% modified by Vincent from Direl's code
ThetaJJ = 2*ThetaJJ - ThetaJJ*Svv*ThetaJJ;  % Unbiased Source Precision Matrix estimator
end
