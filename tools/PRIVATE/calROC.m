function [P, TPR, FPR, TPRL, FPRL, TP, FP, TPL, FPL, TN, FN, TNL, FNL...
    ] = calROC(PMs, PMLaps, GT, nSample, nTH, thMode, sumMode)

[nData, nStat] = size(PMs);
if ~isempty(PMLaps); calLap = 1; [nDataL, ~] = size(PMLaps); else; calLap =0; end
STEP = 1/(nTH-1); P = (1:1:nTH-1)*STEP; P = [0 P];
TPR  = zeros(nTH, nData, nStat); FPR  = zeros(nTH, nData, nStat); 
TP   = zeros(nTH, nData, nStat); FP   = zeros(nTH, nData, nStat); 
TN   = zeros(nTH, nData, nStat); FN   = zeros(nTH, nData, nStat); 
if calLap
    TPRL = zeros(nTH, nDataL);  FPRL = zeros(nTH, nDataL);
    TPL  = zeros(nTH, nDataL);   FPL = zeros(nTH, nDataL);
    TNL  = zeros(nTH, nDataL);   FNL = zeros(nTH, nDataL);
else
    TPRL = [];  FPRL = []; TPL  = [];   FPL = []; TNL  = [];   FNL = [];
end
for i=1:nTH
    
    if calLap
        [PMth, PMLapth] = thresPMall(PMs, PMLaps, nSample, P(i), thMode);
    else
        [PMth, ~] = thresPMall(PMs, {}, nSample, P(i), thMode);
    end
    [ TP(i,:,:), TN(i,:,:), FP(i,:,:), FN(i,:,:)] = ...
        sumStatisticsPMs(PMth, GT, sumMode);
    TPR(i,:,:) = round(TP(i,:,:)./(TP(i,:,:)+FN(i,:,:)), 4);
    FPR(i,:,:) = round(FP(i,:,:)./(FP(i,:,:)+TN(i,:,:)), 4);
    if calLap
        [ TPL(i,:), TNL(i,:), FPL(i,:), FNL(i,:)]     = ...
            sumStatisticsPMs(PMLapth, GT, sumMode);
        TPRL(i,:)  = round(TPL(i,:)./(TPL(i,:)+FNL(i,:)), 4);
        FPRL(i,:)  = round(FPL(i,:)./(FPL(i,:)+TNL(i,:)), 4);
    end
end
end