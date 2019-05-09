function [EEG_IM_awake, ECoG_IM_awake, EEG_IM_ane, ECoG_IM_ane] = formatPMs(EECoGdata, ECoGLap)
% format the PMs for comparison.
% By Vincent 2019.1.9
EEG_IM_awake = EECoGdata(1,1:3); EEG_IM_awake(2,1:2) = EECoGdata(1,4:5); 
EEG_IM_awake(2,3)=ECoGLap(1);
ECoG_IM_awake = EECoGdata(2,1:3); ECoG_IM_awake(2,1:2) = EECoGdata(2,4:5); 
ECoG_IM_awake(2,3)=ECoGLap(1);

EEG_IM_ane = EECoGdata(3,1:3); EEG_IM_ane(2,1:2) = EECoGdata(3,4:5); 
EEG_IM_ane(2,3)=ECoGLap(2);
ECoG_IM_ane = EECoGdata(4,1:3); ECoG_IM_ane(2,1:2) = EECoGdata(4,4:5); 
ECoG_IM_ane(2,3)=ECoGLap(2);
end