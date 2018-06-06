function SIMBIO_make4monkey(Silicone,Inskull,Outskull,Scalp,pos,venant,addpos,Electrodes,addelec,deepelec,elec_labels,id,idelec,outfolder,ratio,v,conds,ow,nobndsplit,scale)
% comment by Vincent
% Silicone, Inskull, Outskull, Scalp: the pathes of the compartments (with vertices and faces)
% pos: Cortex.vertices partion can be 1 or 2
% venant: Venant is true, sources will move to their closest venant nodes
%         EEG 0, and ECoG 1;
% addpos: add position flag, 1 for EEG and ECoG;
% Electrodes: electrodes postion number*3;
% addelec:weather to add electrodes; 
% deepelec: 0 for EEG and 1 for ECoG
% elec_labels: electrodes labels
% id: file name identifier like "Monkey-EEG1"
% idelec: '' for most cases
% outfolder: the output folder "4sh"
% ratio: 1.414, input parameter for tetgen
% v: maximum volume number for finite element method, normally, it is [2 2 2 .1]
% conds: % conductivity for the compartments, normally [.33 .0041 .33 .000001]
% ow:over write flag;
% nobndsplit: ???
% scale: ???

% defalt scale = 1
if nargin < 20
    scale = 1;
end

% configure tetgen path
if isunix
%     tetgen = '/apps/tetgen1.5.0/tetgen';
    tetgen = '/root/Desktop/Datos_Macaque/M.m/external_tool/tetgen1.5.0/tetgen';
%    tetgen = '/root/Desktop/Datos_Macaque/M.m/external_tool/tetgen1.4.3/tetgen';
else
    tetgen = which('tetgen.exe');
end

% used in computing terahedra model
if venant
    warning('if venant is true, sources will move to their closest venant nodes!'); %#ok<WNTAG>
end

% now it is '', and idelec is ''
if addelec
    idelec_mesh = idelec;
else
    idelec_mesh = '';
end

% configure the output file
nmesh = fullfile(outfolder,[id idelec_mesh '_MESHmodel']);
% what is this nmesh?
smesh = fullfile(outfolder,[id idelec_mesh '_MESHmodel.smesh']);
% what is this smesh?
nodefile = fullfile(outfolder,[id idelec_mesh '_MESHmodel.1.node']);
% what is the node file?
elefile = fullfile(outfolder,[id idelec_mesh '_MESHmodel.1.ele']);
% what is the elefile?
facefile = fullfile(outfolder,[id idelec_mesh '_MESHmodel.1.face']);
% what is face file?

if ow || ~(exist(smesh,'file') && exist(nodefile,'file') && exist(elefile,'file') && exist(facefile,'file'))
    % if overwrite flag is 1 or any of the file is missing, recalculate; 
    % surface is the data structure of the headmodel;
    surfs(1) = struct('vertices',Scalp.vertices,'faces',Scalp.faces);
    % surfs(1) scalp + electrodes
    % surfs(2) outskull
    % surfs(3) inskull,            
    % surfs(4+) silicone maybe more than 1ayers)
    if addelec
        surfs(1) = addelectrodes(surfs(1),Electrodes);
        % add electrodes to the surface
    end
    surfs(2) = struct('vertices',Outskull.vertices,'faces',Outskull.faces);
    surfs(3) = struct('vertices',Inskull.vertices,'faces',Inskull.faces);
    facecell = finddisconnsurf(Silicone.faces); % deal with the disconnect surface
    Nfs = numel(facecell); % 4
    for i = 1:Nfs
        surfs(3+i) = struct('vertices',Silicone.vertices,'faces',facecell{i});
        surfs(3+i) = remove_not_refpoints(surfs(3+i)); 
    end
    
    v(5:3+Nfs) = v(4);
    % maximum volume number of silicone for finite element method;
    names = {'scalp','skull','brain'};
    intp = zeros(3+Nfs,3);
    for i = 1:Nfs
        intp(3+i,:) = intpoint(surfs(3+i));
        names(3+i) = {['silicone-' num2str(i)]};
    end
    % silicone names
    obj2 = boundary(patch('vertices',surfs(2).vertices,'faces',surfs(2).faces));
    % obj2 is the boundary of out skull
    obj1 = boundary(patch('vertices',surfs(1).vertices,'faces',surfs(1).faces));
    % obj1 is the boundary of scalp
    [~,indr] = max(pos(:,1));
    intp(3,:) = pos(indr,:);
    [~,indr] = max(surfs(3).vertices(:,1));
    c2 = surfs(3).vertices(indr,:);
    p2 = intersect(obj2,c2,[1e5 0 0],'-all'); 
    p1 = intersect(obj1,c2,[1e5 0 0],'-all');
    close all;
    intp(2,:) = (c2       +p2(1,:))/2;
    intp(1,:) = (p2(end,:)+p1(1,:))/2;
    fsmesh = create_smesh(surfs,names,intp,v,nmesh);
    % creat all the meshes, fmesh
    
    if nobndsplit
        % Y controls if not split the boundary, parameter for tetgen 
        Y = 'Y';
    else
        Y = '';
    end
    
    if addpos
        % generate the .a.node file
        [pp,nn,~] = fileparts(fsmesh);
        posfile = fullfile(pp,[nn '.a.node']);
%         if isunix
%             posfile = fullfile(pp,[nn '.a.node']);
%         elseif ispc
%             posfile = fullfile(pp,[nn '.a.node']);
%         end
        create_node(posfile,pos,[],[]);
        % create node file
        % call tetgen to creat the headmodel, and generate all the files
        % for calculating the leadfield
        sprintf('%s -paACVO7q%fi%s "%s"',tetgen,ratio,Y,fsmesh)
        system(sprintf('%s -paACVO7q%fi%s "%s"',tetgen,ratio,Y,fsmesh));
    else
        system(sprintf('%s -paACVO7q%f%s "%s"',tetgen,ratio,Y,fsmesh));
    end
    % read the results from tetgen output and create the new file.
    [elements,labels] = read_ele(elefile);
    nodes = read_node(nodefile);
    vol = elemvolume(nodes,elements,'signed');
    ind = vol<0; elements(ind,:) = elements(ind, [1 2 4 3]);
    % Correct the 4 ECoG compartments label names into one with label 4.
    labels(labels>4) = 4;
    create_ele(elefile,elements,labels);
    
end

elcfile  = fullfile(outfolder,[id idelec '.elc']);
dipfile  = fullfile(outfolder,[id idelec_mesh '.dip']);
parfile  = fullfile(outfolder,[id idelec_mesh '.par']);
meshfile = fullfile(outfolder,[id idelec_mesh '.v']);

%prepare sensors
if ow || ~(exist(dipfile,'file') && exist(parfile,'file'))
    
    %read nodes and elements
    nodes = read_node(nodefile);
    if ~exist('elements','var')
        [elements,labels] = read_ele(elefile);
    end
    
    disp('writing the dipoles file on disk...')
    if venant
        pos = shift_to_Venant(pos,nodes,elements,labels,3);
    end
    sb_write_dip(pos*scale,dipfile);
    
    % write the parameters file, contains tissues conductivities...
    % grid and simbio call details, mixed together
    disp('writing the parameters file on disk...')
    sb_write_par(parfile,'cond',conds,'labels',unique(labels));
end
if ~exist(meshfile,'file')
    % write the vol.wf in a vista format .v file
    if ~exist('nodes','var')
        nodes = read_node(nodefile);
    end
    if ~exist('elements','var')
        [elements,labels] = read_ele(elefile);
    end

    disp('writing the mesh file on disk...')
    write_vista_mesh(meshfile,nodes*scale,elements,labels);
    
end

disp('writing the electrodes file on disk...')
if ~deepelec
    sb_write_elc(Electrodes*scale,elec_labels,elcfile);
else
    sb_write_elc(Electrodes*scale,elec_labels,elcfile,1);
end

function [pos_venant,inds] = shift_to_Venant(pos,nodes,elements,labels,label)

try %#ok<TRYNC>
    [elements,labels] = read_ele(elements);
    nodes = read_node(nodes);
end

elements  = [elements(:,1),elements(:,2),elements(:,4),elements(:,3)];
label_elements = elements((labels == label),:);
label_nodes_ind = unique(label_elements);

[label_facets, ~] = extract_surface(label_elements);
clear pial_elements
label_surf_nodes_ind = unique(label_facets);

ind_venant_nodes = setdiff(label_nodes_ind,label_surf_nodes_ind);
venant_nodes = nodes(ind_venant_nodes,:);

KDT = KDTreeSearcher(venant_nodes);
near_venant_node_ind = knnsearch(KDT,pos);
pos_venant = venant_nodes(near_venant_node_ind,:);
inds = ind_venant_nodes(near_venant_node_ind);

shift_dist = sqrt(sum((pos_venant-pos).^2,2));

disp('[max, min, median, mean] of distance by which source space nodes where shifted to venant nodes')
disp([max(shift_dist), min(shift_dist), median(shift_dist), mean(shift_dist)]);

function p = intpoint(S)

c = mean(S.vertices);
d = pdist2(c,S.vertices);
[~,ind] = min(d);
o = S.vertices(ind,:);
normals = getNormals(S);
n = normals(ind,:);
obj = boundary(patch('vertices',S.vertices,'faces',S.faces));
p1 = intersect(obj,o, n+eps,'-last'); %#ok<*GTARG>
p2 = intersect(obj,o,-n+eps,'-last'); %#ok<*GTARG>
close all;
p = mean(unique([o;p1;p2],'rows'));

function [S2,electrodes] = addelectrodes(S,electrodes0)

th = .5;
c = mean(S.vertices);
M = [eye(4,3) [-c(:); 1]];
S.vertices = applyM(S.vertices,M);
electrodes0 = applyM(electrodes0,M);
[electrodes,fc] = PontoSurf(electrodes0,S);
indnan = any(isnan(electrodes),2);
electrodes(indnan,:) = electrodes(indnan,:);
P = electrodes(~indnan,:);
fc = fc(~indnan);
%---
TR = triangulation(S.faces,S.vertices);
% original version TriRep, update by vincent
B = cartesianToBarycentric(TR,fc,P);
% old version cartToBary, update by vincent
for i = 1:size(P,1)
    [bmx,indmx] = max(B(i,:));
    if bmx>th
        S.vertices(S.faces(fc(i),indmx),:) = P(i,:);
    end
end
p = unique([electrodes; S.vertices],'rows');
r = sqrt(dot(p,p,2));
S2.vertices = applyM(p,inv(M));
S2.faces = convhulln(p./r(:,[1 1 1]));