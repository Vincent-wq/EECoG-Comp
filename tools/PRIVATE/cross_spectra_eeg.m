function [Svv, Cvv, F] = cross_spectra_eeg(Data, Fs, FMax, nw)
% cross_spectra_eeg use the multitaper methods to estimate the crosspectrum
% of multichannel time series matrix.
%
% Outputs: Svv  : Cross-spectrum matrix [NSensor*NSensor*Nf];
%          Cvv  : Normalized Cross-spectrum matrix [NSensor*NSensor*Nf];
%           F   : The frequency ruler.
% Inputs:  Data : in format of Number of Sensor * Number of time points *
%                 Number of segments;
%           Fs  : Sampling frequency of Data;
%          FMax : The Maximun frequency of the outputs;
%           nw  : time_halfbandwidth, half number of Slepian sequences for 
%                 for the FFT when using multitaper method.
% 
% Created by Pedro Valdes on xxx.xxx.xxx.
% Last modified by Vincent on 27th Oct. 2018

[NSensor, Nt, Nseg] = size(Data); % Get dimensions of the Data .          
tapers  = dpss(Nt,nw); tapers = reshape(tapers,[1,Nt,2*nw]); % generate tappers.
deltaf = Fs/(Nt-1);
F      = (1:(Nt/2))*deltaf;
F      = F(F<FMax);               % The output frequency ruler.
Nf     = length(F);               % Number of frequency points.

Svv = zeros(NSensor,NSensor,Nf);  % init of the cross-spectrum matrix.
Cvv = zeros(NSensor,NSensor,Nf);  % init of the coherence matrix.

for seg = 1:Nseg
    tmp = Data(:, :, seg);
    tmp = repmat(tmp, [1,1,2*nw]).*repmat(tapers, [NSensor,1,1]); 
    % extend the length to improve low frequency estimate
    tmp = fft(tmp,[],2);
    tmp = tmp(:,1:(Nf+1),:); 
    
    for f = 1:Nf
        tmpF=squeeze(tmp(:,f,:)).';
        % cross spectrum calculation for all frequncies
        Svv(:,:,f) = Svv(:,:,f)+cov(tmpF,1); 
    end        
end
Svv   = Svv/Nseg;                % average over segments
for k = 1:Nf
    tmp =squeeze(Svv(:,:,k));
    tmpD = diag(tmp);
    tmpDD = diag(sqrt(tmpD));
    Cvv(:,:,k) = tmpDD\tmp/tmpDD;
end
end