function [Xline,Yline,Zline] = cut_surf(plane,fv)

% This function cuts a surface (defined by vertices and faces as represented 
% by MATLAB patches) by a plane, leading to a curve in 3D space. The
% resulting curve is represented by a set of contigous lines in the space
% 
% Syntax:
% [Xline,Yline,Zline] = cut_surf(plane,SurfData)
% 
% INPUTS:
% plane : A 4-length vector with the parameters of the plane. If plane = [A
% B C D] then every 3D point P = (x,y,z) belonging to the plane satisfies
% plane*[P; 1]' = A*x + B*y + C*z + D = 0
% SurfData : surface structure as represented in MATLAB by patches:
%            SurfData.vertices
%            SurfData.faces
%
% OUTPUTS:
% Xline,Yline,Zline: Matrices with the line coordinates.
% The entire curve can be plotted by simply typing: 
% line(Xline,Yline,Zline,'Properties1',Value1,...);
%
% Pedro Antonio Valdés Hernández
% 
% October 29 2008

warning off %#ok
Xline = []; Yline = []; Zline = [];
oo = ones(size(fv.vertices,1),1);
maxdist = sqrt(max(sum((fv.vertices(fv.faces(:,1),:)-fv.vertices(fv.faces(:,2),:)).^2,2)));
maxdist = max(maxdist,sqrt(max(sum((fv.vertices(fv.faces(:,1),:)-fv.vertices(fv.faces(:,3),:)).^2,2))));
maxdist = max(maxdist,sqrt(max(sum((fv.vertices(fv.faces(:,2),:)-fv.vertices(fv.faces(:,3),:)).^2,2))));
vertx = find(abs(dot([fv.vertices oo],plane(oo,:),2)/norm(plane(1:3)))<maxdist);
indf = ismember(fv.faces,vertx);
[rindf,cindf] = find(indf); %#ok
rindf = unique(rindf);
Nf = length(rindf);
% h = waitbar(0,'cutting surface...');
for i = 1:Nf
   verts = fv.vertices(fv.faces(rindf(i),:),:);
   verts(:,4) = 1;
   diffv(1,:) = diff(verts([1 2],:));
   diffv(2,:) = diff(verts([2 3],:));
   diffv(3,:) = diff(verts([3 1],:));
   alpha = -verts*plane'./(diffv*plane');
   % NaN   : contains
   % < 0   : not contains down
   % -Inf  : parallel down
   % > 1   : not contains up
   % +Inf  : parallel up
   ind = find((alpha<1 & alpha >=0) | (alpha<=1 & alpha >0))  ;
   if ~isempty(ind) && length(ind)==2
       points = verts(ind,1:3) + alpha(ind,[1 1 1]).*diffv(ind,1:3);
       Xline = [Xline points(:,1)]; %#ok
       Yline = [Yline points(:,2)]; %#ok
       Zline = [Zline points(:,3)]; %#ok
   end
%    waitbar(i/Nf,h);
end
% close(h);
warning on %#ok