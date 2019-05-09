function [L,dm] = genlplc(varargin)
%generate the graph laplacian for a spherical head model
% INPUT:
% nneigh: Number of max neighbours for a vertice, default is 6
% opt    : ('standard') [default]  neighbor relations
%          ('dist') Weighted by distance
%           
% Output:
%        L --- laplacian matrix
%       dm --- distance  matrix
% Esin, Andy Hu, 27.03.2016
  
nneigh = 6; opt = 'standard'; % default values 
if     nargin == 1, nneigh = varargin{1};
elseif nargin == 2, nneigh = varargin{1}; opt = varargin{2};
end

sl = importdata('corti869-3000y.dat'); %source locations
sp = dpsetg(sl);
Nv = size(sp,1);

% Calculate the distance matrix
dm = zeros(Nv);
dmvir = dm; % virtual distance matrix that excludes vertice itself
for i = 1:Nv
    d = sqrt(sum((repmat(sl(i,:),Nv,1) - sl(1:Nv,:)).^2,2));   % distance 
    dm(i,:) = d';
    dmvir(i,:) = d'; dmvir(i,i) = Inf;
end

% Laplacian
vxyz = zeros(Nv,nneigh);
L = zeros(Nv);
for i = 1:Nv
    [~,ix] = sort(dmvir(i,:));
    ix = ix(1:nneigh);
    if strcmp(opt,'standard')
        L(i,ix) = -1; L(i,i)  = nneigh;
    else
       vxyz(i,1:nneigh) = ix;
       % L(i,ix) = -sum(DMvir(i,ix))./DMvir(i,ix);
       % L(i,i)  = sum(-L(i,ix));
    end
end

%regular rectangular or circular grids (CNEURO)
if strcmp(opt,'dist')
    L = lplc(vxyz,sl(1:Nv,:));
    L = L (1:3:Nv*3,1:3:Nv*3); % fixed orientation
end

L = sparse(L) + speye(Nv); %make the laplacian full rank
L = L\eye(Nv);
end