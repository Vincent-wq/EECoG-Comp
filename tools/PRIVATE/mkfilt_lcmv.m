function [A]=mkfilt_lcmv(L,C,alpha)
% [A, A1, po]=mkfilt_lcmv(L,C,alpha)
% Calculates the power over all grid points using lcmv
% This is a vector beamformer. For each voxel the dipole direction is chosen, which maximizes power.  
%
% Input:
% L: Lead field matrix, NxMx3 matrix for N channels, M voxels and 3 dipole
%    directions
% C:  NxN matrix for N channels covariance matrix or real part of cross-spectrum
%    (The program uses only the real of C)
% alpha: regularization parameter. In the algorthm C+alpha*eye(N) is
%        inverted. A reasonable value is alpha=trace(C)/(N*100)
%
% Output
% A  : 3-dimension filter
% A1 : 1-dimensional filter along direction with strongest power
% po : Mx1 vector for M voxels, po(i) is the power at the i.th voxel along
%      stronges direction
% 
% History
% by Guido Nolte 
% A1, po removed by Vincent
% License
%   Copyright (C) 2015
%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with this program.  If not, see http://www.gnu.org/licenses/.

C=real(C);

if nargin<3
    alpha=.05*trace(C)/length(C);
end
[nchan, ns, ndim]=size(L);
Cr=C+alpha*eye(nchan);

Crinv=inv(Cr);

A=zeros(nchan,ns,ndim);
for i=1:ns
    Lloc=squeeze(L(:,i,:));
    A(:,i,:)=reshape((inv((Lloc'*Crinv*Lloc))*Lloc'*Crinv)',nchan,ndim); %#ok<*MINV>
end
end
%% removed by Vincent 
% po=zeros(ns,1);
% A1=zeros(nchan,ns);
% for i=1:ns
%     Aloc=transpose(squeeze(A(:,i,:)));
%     Ploc=Aloc*C*Aloc';
%     [u,s,~]=svd(Ploc);
%     po(i)=s(1,1);
%     A1(:,i)=Aloc'*u(:,1);
% end

