function map = Bipolar3(nbins)
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
if (nbins < 13)
    error('The number of bins should be >= 13!');
end
if (abs(nbins - 2*floor(nbins/2)) ~= 1)
    error('Enter an odd number!');
end
n = ceil(nbins/2);
bipolar = [255 255   0; ...
           255 170  43; ...
           255  85  85; ...
           213  43  43; ...
           150   0   0; ...
           160 128 128; ...
           170 170 170];
xx = linspace(1,7,n);
map = interp1(1:7, bipolar, xx, 'linear');
map = [map(1:end-1,:); rot90(map,2)];
map = flipud(map)/255;