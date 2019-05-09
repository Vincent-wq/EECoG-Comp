function [data, ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglabReader(NAME,ALLEEG,EEG,CURRENTSET,ALLCOM)
EEG = pop_loadset( 'filename',[NAME '.set']);
[ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
data = EEG.data;
eeglab redraw
end