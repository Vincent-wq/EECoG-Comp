function [segData] = def_segment(Data, Fs, deltaF)
% input: Data  : Time series matrix: nSensor * NtTotal;
%        Fs    : Sampling frequency;
%        deltaF: Frequency resolution;
% output data  : Segmented data: nSensor * NtTotal * NSeg
%                by Vincent 2018.4.27
[nSensor, NtTotal] = size(Data);
NpSeg  = round(Fs/deltaF);
NSeg   = floor(NtTotal/NpSeg); Data = Data(:,1:NpSeg*NSeg);
segData = reshape(Data, [nSensor, NpSeg, NSeg]);
end