function [ALLEEG, EEG, CURRENTSET] = eegLaplacianFormat(EEG, ALLEEG, M, showMAP)
if showMAP
    MapMontage(M)
end
[G,H] = GetGH(M);
X = CSD (EEG.data, G, H);
EEG.data = X;
EEG.setname    = [EEG.setname '-Lap']; %data.cfg.dataset;
EEG.comments   = 'preprocessed with CSD';
EEG.nbchan     = size(X,1);
EEG.trials     = size(1,2);
EEG.pnts       = size(X,2);
[ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG);
end