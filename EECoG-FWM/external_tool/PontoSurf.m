function [pn,faces,normals] = PontoSurf(p,S,c,sink)
% comment by Vincent
% what is this function for?
% Inputs:
% p    : 
% S    : Seems the surface with .vertices and surfaces
% c    : 
% sink : 
%
if nargin > 2
    S.vertices = S.vertices-c(ones(size(S.vertices,1),1),:);
    p = p-c(ones(size(p,1),1),:);
end
if nargin < 4
    sink = 0;
end

f = figure('Visible','off'); h = patch(S,'Visible','off');
obj = boundary(h);
norms = sqrt(dot(p,p,2));
normals = p./norms(:,[1 1 1]);
vector = (2*max(norms))*normals;
[coll,rvalue,fc,nrms] = intersect(obj,zeros(size(p,1),3),vector,'-all');
if ~iscell(coll)
    coll = {coll};
    rvalue = {rvalue};
    fc = {fc};
    nrms = {nrms};
end
close(f);

% this must be fixed
% English explaination?
siz = cellfun('size',coll,1);
pn = zeros(size(p,1),3);
faces = zeros(size(p,1),1);
normals = zeros(size(p,1),3);
for i = 1:size(p,1)
    if siz(i)>1
        d = coll{i}-p(ones(siz(i),1),:);
        d = dot(d,d,2);
        [~,ind] = min(d);
        if sink
            pn(i,:) = rvalue{i}(ind)*(1-sink)*vector(i,:);
        else
            pn(i,:) = coll{i}(ind,:);
        end
        faces(i) = fc{i}(ind);
        normals(i,:) = nrms{i}(ind);
    elseif siz(i)==1
        if sink
            pn(i,:) = rvalue{i}*(1-sink)*vector(i,:);
        else
            pn(i,:) = coll{i};
        end
        faces(i) = fc{i};
        normals(i,:) = nrms{i};
    else
        pn(i,:) = NaN(1,3);
        faces(i) = NaN;
        normals(i,:) = NaN(1,3);
    end
end

if nargin > 2
    pn = pn+c(ones(size(p,1),1),:);
end