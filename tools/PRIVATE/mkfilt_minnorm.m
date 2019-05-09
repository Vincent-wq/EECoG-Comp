function A=mkfilt_minnorm(L,alpha0)
% makes spatial filter according to minimum norm solution
% usage  A=mkfilt_minnorm(L); or  A=mkfilt_minnorm(L,alpha0);
%
% input L:  NxMxP leadfield tensor for N channels, M voxels, and 
%           P dipole directions. Typically P=3. (If you do MEG for 
%           a spherical volume conductor or reduce the rank, you must 
%           reduce L such that it has full rank for each voxel, such that,
%           e.g., P=2)
%       alpha: optional regularization parameter (default is .05 corresponding 
%             to 5% of the average of the eigenvalues.) 
% 
% output A: NxMxP tensor of spatial filters. If x is the Nx1 data vector at time t. 
%           then A(:,m,p)'*x is the source activity at time t in voxel m in source direction
%           p. 
% modified by Vincent for accuracy of division.
if nargin<2
    alpha0=.05;
end

[nchan,ng,ndum]=size(L);
LL=reshape(L,nchan,ng*ndum);
R=LL*LL';
alpha=alpha0*trace(R)/length(R);
A=reshape((LL'/(R+alpha*eye(nchan)))',nchan,ng,ndum);
return