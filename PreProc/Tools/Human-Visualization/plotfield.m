function plotfield(fig,alpha,cam,linestyle,Cortex,points,field,cl)

% fig: opens figure
% Cortex: patch
% points: ej. patch.vertices(ind,:) or []
% field: vector same length of size(patch.vertices,1)
% cl: guiate por el switch de abajo

if isempty(points)
    points = Cortex.vertices;
end
indp = ismember(Cortex.vertices,points,'rows');
values = zeros(size(Cortex.vertices,1),1);
values(indp) = field;

% [ut,i1t,it2] = unique(values(indp)); %#ok
[ut,i1t,it2] = unique(values);
switch cl
    case 'gray'
        col = gray(size(ut,1));
    case 'hot'
        col = hot(size(ut,1));
    case 'hsv'
        col = hsv(size(ut,1));
    case 'jet'
        col = jet(size(ut,1));
    case 'cool'
        col = cool(size(ut,1));
    case 'bone'
        col = bone(size(ut,1));
    case 'pink'
        col = pink(size(ut,1));
    case 'winter'
        col = winter(size(ut,1));
    case 'autumn'
        col = autumn(size(ut,1));
    case 'summer'
        col = summer(size(ut,1));
    case 'spring'
        col = spring(size(ut,1));
    case 'copper'
        col = copper(size(ut,1));
    case 'spectral'
        col = spectral(size(ut,1));
    case 'varycolor'
        col = varycolor(size(ut,1));
end
FaceVertexCData = zeros(size(Cortex.vertices));
FaceVertexCData0 = col(it2,:);
FaceVertexCData(indp,:) = FaceVertexCData0(indp,:);

Cortex.FaceVertexCData = FaceVertexCData;
plotsurf(Cortex,fig,alpha,cam,linestyle,FaceVertexCData);