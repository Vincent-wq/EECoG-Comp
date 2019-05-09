function plotelec(S,proj,pos,names,radius,height,color,c,propname)

if nargin == 9
   dis = propname{1};
   siz = propname{2};
   wei = propname{3};
end
%---
[X,Y,Z] = cylinder(radius,10);
Z(2,:) = height;
se1 = surf2patch(X,Y,Z);
%---
theta = pi/4;
[X,Y,Z] = sphere(10);
se2 = surf2patch(X,Y,Z);
[~,T,~] = cart2sph(se2.vertices(:,1),se2.vertices(:,2),-se2.vertices(:,3));
inds = find((pi/2-T)>theta);
se2 = opensurf(se2,inds);
r = radius/sin(theta);
vn = se2.vertices./repmat(sqrt(dot(se2.vertices,se2.vertices,2)),1,3);
se2.vertices = vn.*r;
se2.vertices(:,3) = se2.vertices(:,3) + height + radius*cot(theta);
se.vertices = [se1.vertices; se2.vertices];
se.faces = [se1.faces; se2.faces+size(se1.vertices,1)];
%---
[P,T,~] = cart2sph(se2.vertices(:,1),se2.vertices(:,2),se2.vertices(:,3)-(height+radius*cot(theta)));
thetamax = max(T(:));
ind = find(round(T*1000)==round(thetamax*1000));
P = P(ind);
[~,inds] = sort(P);
P = P(inds);
[X,Y,Z] = sph2cart(P,(theta-pi/2)*ones(size(P,1),1),r*ones(size(P,1),1));
bpoints2(:,1) = X;
bpoints2(:,2) = Y;
bpoints2(:,3) = Z+height+radius*cot(theta);
bedges = ind(inds)+size(se1.vertices,1);
bedges = [bedges' bedges(1)];
Nv = size(se.vertices);
bedges2 = [(Nv+1):(Nv+numel(bedges)-1) Nv+1];
for i = 1:numel(bedges)-1
    face = [bedges(i:i+1) bedges2(i+1:-1:i)];
    se.faces = [se.faces; face];
end
se.vertices = [se.vertices; bpoints2];
%---
if proj
    [pos,tri] = project_points_surf(S,pos,c);
else
    [~,tri] = project_points_surf(S,pos,c);
end
sep = se;
for i = 1:size(pos,1)
    p1 = S.vertices(tri(i,1),:);
    p2 = S.vertices(tri(i,2),:);
    p3 = S.vertices(tri(i,3),:);
    n = [p1;p2;p3]\ones(3,1);
    n = n/norm(n);
    R = makeR(n);
    M = [R pos(i,:)'; 0 0 0 1];
    sep.vertices = applyM(se.vertices,M);
    plotsurf(sep,0,1,1,'none',color);
    if ~isempty(names) && nargin == 9
        name = names{i};
        if isnumeric(name)
            name = num2str(name);
        end
        h = text(pos(i,1)+dis*n(1),pos(i,2)+dis*n(2),pos(i,3)+dis*n(3),name);
        set(h,'Color',[0 0 1]);
        set(h,'FontSize',siz);
        set(h,'FontWeight',wei);
        set(h,'HorizontalAlignment','center');
    end
end

function S = opensurf(S,inds)

for i = 1:length(inds)
    [r{i},c] = find(S.faces == inds(i)); %#ok
end
S.faces(unique(cat(1,r{:})),:) = [];
%---
S.vertices(inds,:) = [];

tmpfaces = S.faces;
for i = 1:length(inds)
    ind = find(S.faces >= inds(i));
    tmpfaces(ind) = tmpfaces(ind)-1;
end
S.faces = tmpfaces;

function [Pn,Tn] = project_points_surf(S0,Pn,c)

Ne = size(Pn,1);
St = S0;
St.vertices = S0.vertices-c(ones(size(S0.vertices,1),1),:);
Pnt = Pn - c(ones(Ne,1),:);
Tn = zeros(size(Pn,1),3);
hwait = waitbar(0,'re-projecting electrodes...');
for i = 1:Ne
    p0 = Pnt(i,:);
    d = St.vertices-p0(ones(size(St.vertices,1),1),:);
    d = sqrt(dot(d,d,2));
    [~,ind0] = sort(d);
    g = []; p = []; T = [];
    Ni = 20;
    for k = 1:10
        ind = ind0(1:Ni);
        tri = get_tri(St,ind);
        m = length(tri);
        for j = 1:m
            p1 = St.vertices(tri(j,1),:)';
            p2 = St.vertices(tri(j,2),:)';
            p3 = St.vertices(tri(j,3),:)';
            prms = [p1-p2 p1-p3 p0']\p1;
            if all(0<=prms) && all(prms([1 2])<=1) && sum(prms([1 2]))<=1
                g = [g; prms(3)]; %#ok<AGROW>
                p = [p; prms(3)*p0]; %#ok<AGROW>
                T = [T; tri(j,:)]; %#ok<AGROW>
            end
        end
        if ~isempty(g)
            break;
        else
            Ni = Ni*10;
        end
    end
    [~,ind] = max(g);
    Pn(i,:) = c + p(ind,:);
    Tn(i,:) = T(ind,:);
    waitbar(i/Ne,hwait);
end
close(hwait);