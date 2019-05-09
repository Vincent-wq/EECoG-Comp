function [hFig] = patchVis(P, Nodes, Val, conf)
%       input: 
%                   P: patch object to be plotted on;
%               Nodes: nodes of the patches Ns*3 
%                 Val: data vector to plot, Nsl*3
%                conf: cofiguarations for the visualization.
%                      more specs.: ...
%      output:   hFig: the handle of the figure.
%                   
%  v1.01   By Vincent @ 1st Nov. 2018

%% input process
% match node of input
if isempty(Nodes)
    Nodes = P.Vertices;
end
% format the input values
indp = ismember(P.Vertices,Nodes,'rows');
ValDisp = zeros(size(P.Vertices,1),1)+min(Val);

if isfield(conf,'log') && conf.log ==1
    ValDisp(indp) = log10(Val);
    if ~isfield(conf,'clim') 
        if isfield(conf,'logCut')
            conf.clim=[min(ValDisp)+conf.logCut max(ValDisp)];
        else
            conf.clim=[min(ValDisp) max(ValDisp)];
        end
    else
        conf.clim=[log10(conf.clim(1)) log10(conf.clim(2))];
    end
else
    ValDisp(indp) = Val;
    if ~isfield(conf,'clim') 
        conf.clim=[min(ValDisp) max(ValDisp)];
    end
end
% test
%% patch initialization
if conf.ifFig, figure; else hold on; end
%%
hP = patch('Vertices', P.Vertices, 'Faces', P.Faces, ...
      'FaceVertexCData',ValDisp,...
      'FaceColor', 'interp',...
      'FaceLighting','phong',...
      'LineStyle',conf.linestyle,...
      'LineWidth',conf.LineWidth,...
      'EdgeColor', 'interp',...
      'FaceAlpha',conf.alpha);
axis equal;
shading interp

%% light & lighting
lightangle(-45,30)
hP.FaceLighting = 'gouraud';    
hP.AmbientStrength = conf.material(1);
hP.DiffuseStrength = conf.material(2);
hP.SpecularStrength = conf.material(3);
hP.SpecularExponent = conf.material(4);
hP.SpecularColorReflectance =  conf.material(5);
hP.BackFaceLighting = 'unlit';
if isfield(conf,'FValpha')
    hp.FaceVertexAlphaData  = conf.FValpha;
end
camlight(0,180);
camlight(0,0);
camlight('headlight');
camlight('right');
camlight('left');
title(conf.TITLE);
%%set colorbar
hFig = gca;
colormap(conf.colormap);
if isfield(conf,'ifColorBar') && conf.ifColorBar == 1
    colorbar('peer',hFig);
    set(hFig,'CLim',conf.clim)
end
%% view and tools
view(-151,30); 
cameratoolbar;
end