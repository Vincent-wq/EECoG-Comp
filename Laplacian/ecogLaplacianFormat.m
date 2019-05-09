function [ALLEEG, EEG, CURRENTSET] = ecogLaplacianFormat(EEG, ALLEEG, Lap)
% calculate the ECoG laplacian based on the cortex surface
% By Vincent&Pedro @ 2018.11.12
EEG.data = Lap;
EEG.setname    = [EEG.setname '-Lap']; %data.cfg.dataset;
EEG.comments   = 'preprocessed with cortex Laplacian';
EEG.nbchan     = size(Lap,1);
EEG.trials     = size(1,2);
EEG.pnts       = size(Lap,2);
[ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG);
end