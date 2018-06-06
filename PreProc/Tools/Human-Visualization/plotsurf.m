function plotsurf(fv,iffig,alpha,lighton,linestyle,color,coloredge)

% fv: patch
% iffig: new figure 1
% lighton: light 1
% linestyle: '-' wireframe, 'none' none
% color: [1 0 0]
% coloredge: [1 0 0] o 'r'?

if isempty(linestyle)
    linestyle = 'none';
end

if nargin == 5 || nargin == 6
    coloredge = [.3 .3 .3];
end
% colordef black
if ischar(fv)
    fv = load(deblank(fv));
    fv = fv.SurfData;
end
V = fv.vertices; F = fv.faces;
if nargin >= 6
    if size(color,1)==1,
        colors = color(ones(size(V,1),1),:);
    elseif size(color,1)==size(V,1)
        colors = color;
    end
else
    h = figure; normals = get(patch(fv),'VertexNormals'); close(h)
    colors = abs(normals./(repmat(sqrt(dot(normals,normals,2)),1,3)+eps));
    % colors = abs(V./(repmat(sqrt(dot(V,V,2)),1,3)+eps));
end
if nargin == 4 || isempty(linestyle)
    linestyle = '-';
end
if iffig, figure; else hold on; end
patch('Vertices', V, 'Faces', F, ...
      'FaceVertexCData',colors,...
      'FaceColor', 'interp',...
      'FaceLighting','phong',...
      'LineStyle',linestyle,...
      'LineWidth',1,...
      'EdgeColor',coloredge,...
      'AmbientStrength',.4,...
      'FaceAlpha',alpha);
axis equal
% camlight('right')%,'infinite')
% camlight('left')%,'infinite')
if iffig && lighton,
    camlight(0,180)%,'infinite')
    camlight(0,0)%,'infinite')
end
cameratoolbar