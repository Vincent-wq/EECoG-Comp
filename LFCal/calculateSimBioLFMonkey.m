%% cleaning working space and configure the working directory
clc; clear;
% working path
folder = '/root/Desktop/Datos_Macaque/';
datafolder = fullfile(folder,'newSimBioModel');
% datafolder = '/root/Desktop/Datos_Macaque/newSimBioModel';
meshFolder = fullfile(datafolder,'4sh');
% add all the external toolboxes needed
% addpath(genpath(fullfile(folder,'M.m','external_tool')));
addpath(genpath(fullfile(folder,'M.m','external_tool','iso2mesh')));
% Simbio conf
% Programfolder = '/root/Desktop/Datos_Macaque/M.m/external_tool/simbio';
SimbioFolder  = fullfile(folder, 'M.m', 'external_tool', 'simbio');
LF_outFolder  = fullfile(datafolder, 'LF');
compList      = fullfile(datafolder, 'compList.txt');
idelec    = '';
ow_simbio = true;
sh_only   = false;
verSimbio = 'NoOutput';
%% Set head model and tetrahedralize
ratio = 1.414;
v     = [10 2 5 .1];
conds = [.33 .0041 .33 .000001];
% Tetgen parameters
nobndsplit = false;
ow_eeg  = true;
ow_ecog = true;
%% loading headmodel elements
cortex = fullfile(datafolder, 'Cortex-mid_Su.mat');
silicone = fullfile(datafolder, 'EcoG-elecs_Su-silicone_sheet.mat');
inskull = fullfile(datafolder, 'Inskull_Su-corrected.mat');
outskull = fullfile(datafolder, 'Outskull_Su-corrected.mat');
scalp = fullfile(datafolder, 'Scalp_Su.mat');
EEG_elec_file = fullfile(datafolder, 'EEG-elecs_Su.mat');
EcoG_elec_file = fullfile(datafolder, 'EcoG-elecs_Su-embbed.mat');
% load and set cortex resolution
CortexResolution = 2000;
Cortex = load(cortex);
Cortex_red = reducepatch(Cortex,CortexResolution/size(Cortex.vertices,1));
[Cortex_red.vertices,Cortex_red.faces] = meshcheckrepair(Cortex_red.vertices,Cortex_red.faces,'meshfix');
save( fullfile(datafolder, ['Cortex_red_',num2str(CortexResolution),'_vincent.mat']), 'Cortex_red');
% load silicone
Silicone = load(silicone);
facecell = finddisconnsurf(Silicone.faces);
N = numel(facecell);
Si = struct;
for i = 1:N
    [Si(i).vertices,Si(i).faces] = meshcheckrepair(Silicone.vertices,facecell{i},'meshfix');
end
Silicone = merge_patches(Si);
% load inskull
Inskull = load(inskull);
[Inskull.vertices,Inskull.faces] = meshcheckrepair(Inskull.vertices,Inskull.faces,'meshfix');
% check and fix inskull and cortex
C = Cortex_red;
m = mean(C.vertices);
m = m(ones(size(C.vertices,1),1),:);
C.vertices = C.vertices-m;
C.vertices = C.vertices*1.02;
C.vertices = C.vertices+m;
[S.vertices,S.faces] = surfboolean(Inskull.vertices, Inskull.faces(:,[1 3 2]),'or',C.vertices,C.faces(:,[1 3 2]));
S.vertices = S.vertices(:,1:3);
S.faces = S.faces(:,1:3);
[S.vertices,S.faces] = meshcheckrepair(S.vertices,S.faces,'meshfix');
Inskull = S;
% load outskull
Outskull = load(outskull);
[Outskull.vertices,Outskull.faces] = meshcheckrepair(Outskull.vertices,Outskull.faces,'meshfix');
% load skull
Scalp = load(scalp);
[Scalp.vertices,Scalp.faces] = meshcheckrepair(Scalp.vertices,Scalp.faces,'meshfix');
pos = Cortex_red.vertices;
% load EEG electrodes
load(EEG_elec_file);
eeg_elec = PontoSurf(electrodes,Scalp,mean(Inskull.vertices)); 
clear electrodes;
eeg_elec_labels = create_numeric_labels(size(eeg_elec,1));
% load ECoG electrodes
load(EcoG_elec_file);
N = size(XYZ,1);
ecog_elec_labels = create_numeric_labels(N);
d = pdist2(XYZ,Silicone.vertices);
[~,ind] = sort(d,2);
ind = ind(:,1:2);
P1 = Silicone.vertices(ind(:,1),:);
P2 = Silicone.vertices(ind(:,2),:);
C = mean(Inskull.vertices);
C = C(ones(N,1),:);
V1 = P1-C;
V2 = P2-C;
D1 = dot(V1,V1,2);
D2 = dot(V2,V2,2);
for i = 1:N
    if D1(i)>D2(i)
        XYZ(i,:) = P2(i,:); 
    else 
        XYZ(i,:) = P1(i,:); 
    end
end
ecog_ele = XYZ; clear XYZ;
%% head modeling
% EEG head model
SIMBIO_make4monkey(Silicone,Inskull,Outskull,Scalp,pos,0,1,eeg_elec,1,0,eeg_elec_labels,...
    'Monkey-EEG','',meshFolder,ratio,v,conds,ow_eeg,nobndsplit);
% ECoG head model
SIMBIO_make4monkey(Silicone,Inskull,Outskull,Scalp,pos,0,1,ecog_ele,0,1,ecog_elec_labels,...
    'Monkey-ECoG','',meshFolder,ratio,v,conds,ow_ecog,nobndsplit);
% Simbio to compute leadfield
SIMBIO_for(meshFolder, SimbioFolder, LF_outFolder, compList,idelec, ow_simbio, sh_only, verSimbio)


