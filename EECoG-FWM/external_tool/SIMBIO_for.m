function SIMBIO_for(DATAFolder,Programfolder,OUTFolder,list,idelec,ow,sh_only,ver)
% comment by Vincent
% DATAFolder: Input file path
% Programfolder: simbio path 
% OUTFolder: leadfield output folder
% list: file of variable name list
% idelec: electrode id
% ow: over write flag
% sh_only: ???
% ver: simbio version 'NoOutput'

%% read name list into ids
tline = 1;
c = 1;
fid = fopen(list);
while tline ~= -1
    tline = fgetl(fid);
    ids{c} = tline; 
    c = c+1;
end
fclose(fid);
%% 
ids = char(ids{1:end-1});
N = size(ids,1);
if nargin == 8
    switch ver
        % simbio version selection
        case 'Venant'
            Program = fullfile(Programfolder,'ipm_linux_opt_Venant');
        case 'NoOutput'
            Program = fullfile(Programfolder,'ipm_linux_opt_NoOutput');
    end
else
    % default simbio version
    Program = fullfile(Programfolder,'ipm_linux_opt_Venant');
end

%% calculate all the leadfields
for i = 1:N
    ID = deblank(ids(i,:));
    % input folder and data
    meshfile = fullfile(DATAFolder,[ID '.v']);
    elcfile = fullfile(DATAFolder,[ID idelec '.elc']);
    dipfile = fullfile(DATAFolder,[ID '.dip']);
    % reconsfile = fullfile(DATAFolder,[ID '.dip.mat']);
    outfile = fullfile(OUTFolder,[ID idelec]);
    % matfile: the output file
    matfile = fullfile(OUTFolder,[ID idelec '.mat']);
    parfile = fullfile(DATAFolder,[ID '.par']);
    logfile = fullfile(OUTFolder,[ID idelec '.log']);
    % the bash file
    shfile =  fullfile(OUTFolder,[ID idelec '.sh']);
    asciifile = fullfile(DATAFolder,[ID '.v.ascii']);
    potfile = fullfile(DATAFolder,[ID '.v.potential']);
    %% generate the bash file
    cmd = sprintf('%s -i sourcesimulation -h %s -s %s -dip %s -o %s -p %s -fwd FEM -sens EEG 2>&1 > %s',...
        Program,meshfile,elcfile,dipfile,outfile,parfile,logfile);
    fid = fopen(shfile,'wt');
    fprintf(fid,'#!/usr/bin/env bash\n\n');
    fprintf(fid,cmd);
    fclose(fid);
    
    %% execute the command and calculate leadfield
    if (ow && ~sh_only) || ~(exist(outfile,'file') || sh_only)
        [s,r] = unix(['sh ' shfile]);
        disp('simbio out:')
        disp(s)
        disp(r)
        %delete old file if recalculate
        if exist(matfile,'file')
            delete(matfile)
        end
    end
    if ~exist(matfile,'file') && ~sh_only
        disp(outfile)
        K = sb_read_msr(outfile);
%   what is the reconsfile for ?? vincent
%         if exist(reconsfile,'file')
%             load(reconsfile);
%             dip3rec = 3*(dip_recons-1);
%             dip_recons1 = dip3rec+1;
%             dip_recons2 = dip3rec+2;
%             dip_recons3 = dip3rec+3;
%             dip_recons123 = [dip_recons1(:) dip_recons2(:) dip_recons3(:)]';
%             dip_recons123 = dip_recons123(:)';
%             K = K(:,dip_recons123); %#ok<NASGU>
%         end
        saveK(matfile,K)
    end
    if exist(asciifile,'file')
        delete(asciifile);
    end
    if exist(potfile,'file')
        delete(potfile);
    end
end

function saveK(matfile,K)

save(matfile,'K');
