function [svName]=ecogGeo4Laplacian(eleFile,cortexFile,svName)
% calculate the ECoG laplacian based on the cortex surface
% By Vincent&Pedro @ 2018.11.12

% for test
% eleFile = 'ECoG_ele128'; cortexFile = 'Cortex_red_10007';
%% load geometry of cortex and electrodes
% load and remove bad ecog electrodes
% ECoG_ele = load(eleFile, 'ECoG_ele');                                  % load electrode profile
rmInd = [38 50 93];                                   % bad electrdos to remove
ECoG_ele.pos = ECoG_ele.ECoG_ele.pos; ECoG_ele.pos(rmInd,:)=[]; 
ECoG_ele.label = ECoG_ele.ECoG_ele.Labels;  ECoG_ele.label(rmInd)=[];
ECoG_ele=rmfield(ECoG_ele,'ECoG_ele');
[nChan,~]=size(ECoG_ele.pos);
% load format cortex patch
cortex = load(cortexFile,'Cortex_red');                 % load cortex profile
cortex.Vertices = cortex.Cortex_red.vertices;
cortex.Faces = cortex.Cortex_red.faces;
cortex = rmfield(cortex, {'Cortex_red'});
[nS,~]   =size(cortex.Vertices);

%%  locate the left hemisphere and get the inds, check with vis
inds = find(cortex.Vertices(:,1)<mean(cortex.Vertices(:,1)));
conf = struct('colormap',{'parula'},...
    'ifFig',{1},...
    'linestyle',{'none'},...
    'LineWidth',{4},...
    'alpha',{0.7},...
    'TITLE',{'test Image'},...
    'log',{1},...
    'material',{[0.6, 0.5, 0.3, 1, 0.6]});
patchVis(cortex, cortex.Vertices(inds,:), ones(length(inds),1), conf);
nsL = length(inds);             % number of sourcess on the left cortex
dchno=pdist2(ECoG_ele.pos, cortex.Vertices(inds, :), 'seuclidean'); 
[m,j]=min(dchno,[],2);
[C,ia,ic] = unique(j);
S = sparse(nChan, nsL);
for i = linspace(1,nChan,nChan)
    S(i,j(i))=1;
end
% calculate cortex laplacian
[Lo, ~] = mesh_laplacian(cortex.Vertices,cortex.Faces);
L = Lo(inds,inds);          % Take left hemisphere
% savae results
if ischar(svName)
    save(svName,'S','L','nChan','nS','nsL','inds','j','cortex','ECoG_ele');
end
end 
% check for lefthemisphere difference
% k=0;
% for i = 1: length(inds)
%     if abs(sum(Lo(inds(i),:))-sum(L(i,:)))>0.0001
%         disp(['ind: ',num2str(inds(i)),', Lo-L: ',num2str(abs( sum(Lo(inds(i),:))-sum(L(i,:)) ))])
%         k=k+1;
%         hold on
%         scatter3(cortex.Vertices(inds(i),1), cortex.Vertices(inds(i),2), cortex.Vertices(inds(i),3),'r')
%     end
% end
% k