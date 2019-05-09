function plottetra(fv,iffig,alpha,lighton,linestyle,color,coloredge)

% fv: patch
% iffig: new figure 1
% lighton: light 1
% linestyle: '-' wireframe, 'none' none
% color: [1 0 0]
% coloredge: [1 0 0] o 'r'?

if nargin == 5 || nargin == 6
    coloredge = [.3 .3 .3];
end

fv1 = struct('vertices',fv.nodes,'faces',fv.elements(:,[1 2 3]));
fv2 = struct('vertices',fv.nodes,'faces',fv.elements(:,[1 3 4]));
fv3 = struct('vertices',fv.nodes,'faces',fv.elements(:,[2 3 4]));
fv4 = struct('vertices',fv.nodes,'faces',fv.elements(:,[1 2 4]));

plotsurf(fv1,iffig,alpha,lighton,linestyle,color,coloredge);
plotsurf(fv2,0,alpha,0,linestyle,color,coloredge);
plotsurf(fv3,0,alpha,0,linestyle,color,coloredge);
plotsurf(fv4,0,alpha,0,linestyle,color,coloredge);