function map = bipolar2(nbins)
% map = colormapBipolar;
% This function creates a bipolar colormap
%
% INPUT:
% nbins:    The number of bins for colormap; by default nbins = 97. This is
%           restricted to be an odd number >= 13.
%
% OUTPUT:
% map:      Bipolar colormap; dimension is nbins x 3

if (nargin == 0), nbins = 97; end
% if (nbins < 13)
%     error('The number of bins should be >= 13!');
% end
% if (abs(nbins - 2*floor(nbins/2)) ~= 1)
%     error('Enter an odd number!');
% end
n = ceil(nbins/2);
bipolar = [255 255   0; ...
           255 170  43; ...
           255  85  85; ...
           230 43  43; ...
           200   0   0; ...
           230 128 128; ...
           255 255 255];
       
xx = linspace(1,7,n);
map = interp1(1:7, bipolar, xx, 'linear');
map = map ./ 255;
map = flipud(map);

cm = [   0         0    1.0000
         0    0.1000    0.9000
         0    0.2000    0.8000
    0.3000    0.4000    0.9000
    0.1000    0.4000    0.8000
    0.2000    0.6000    0.8000
    0.7000    0.8000    0.9000];
    
neg = interp1(1:size(cm,1), cm, linspace(1,size(cm,1),n-1), 'linear');
map = [neg; map];


% map = [map(1:end-1,:); rot90(map,2)];
% map = flipud(map)/255;


% bipolar = [255 255   0; ...
%            255 170  43; ...
%            255  85  85; ...
%            213  43  43; ...
%            170   0   0; ...
%            213 128 128; ...
%            255 255 255];
