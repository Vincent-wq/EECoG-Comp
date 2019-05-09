function show_surf_over_image(image,surfs,slices,Npoints,color)

% This function shows a set of surfaces over an image volume (NIFT or
% ANALYZE). This function needs SPM 
% SYNTAX:
% show_surf_over_image(image,surfs,slices)
% image: pathname to the NIFTI file
% surfs: column vector of pathnames to the surfaces, as built using surfs = char(surf1,surf2,...)
% slices: array specifying the slices in the NIFTI file, for example for sagittal slices:
%         slices = [s1 0 0;...
%                   s2 0 0;
%                   ...
%                   s3 0 0];
%         for coronal slices;
%         slices = [0 s1 0;...
%                   0 s2 0;
%                   ...
%                   0 s3 0];
%         for axial slices;
%         slices = [0 0 s1;...
%                   0 0 s2;
%                   ...
%                   0 0 s3];       

warning('off','MATLAB:conversionToLogical');
ns = size(slices,1);
nc = ceil(sqrt(ns));
nr = ceil(ns/nc);
V = spm_vol(image);
I = spm_read_vols(V);
I(I>prctile(I(:),99)) = prctile(I(:),99);
I = permute(I,[2 1 3]); %I = zeros(size(I));
v = 90*eye(3);
if ~exist('color','var')
    color = 'yrgbcmyrgbcm';
end
figure
h = waitbar(0,'printing overlaid surface...');
for i = 1:ns 
    hold on
    if ns > 1
        subplot(nr,nc,i);
    end
    slicei = slices(i,:);
    fig = slice(I,slicei(1),slicei(2),slicei(3));
    slicei = logical(slicei);
    fign = fig(~slicei(1:3));
    figy = fig(slicei(1:3));
    set(fign(1),'Visible','off');
    set(fign(2),'Visible','off');
    set(figy,'LineStyle','none');
    colormap gray; axis equal; hold on; view(v(slicei,[1 3])); axis off
    plane = double(slicei);
    plane(4) = -sum(slices(i,:));
    if exist('Npoints','var') && numel(Npoints)==1
        Npoints = Npoints(ones(1,size(surfs,1)));
    end
    for s = 1:size(surfs,1)
        if ischar(surfs)
            Surf = load(deblank(surfs(s,:)));
        elseif isstruct(surfs)
            Surf = surfs(s);
        end            
        if isfield(Surf,'Surf')
            Surf = Surf.Surf; 
        end
        if ~isfield(Surf,'SurfData')
            Surf.SurfData = Surf;
            Surf = rmfield(Surf,'vertices');
            Surf = rmfield(Surf,'faces');
        end
        if exist('Npoints','var') && Npoints(s)
            Surf.surfData = reducepatch(Surf.surfData,Npoints(s)/size(Surf.SurfData.vertices,1));
        end
        Surf.SurfData.vertices = [Surf.SurfData.vertices ones(size(Surf.SurfData.vertices,1),1)]*inv(V.mat)';
        Surf.SurfData.vertices(:,4) = [];
        [Xline,Yline,Zline] = cut_surf(plane,Surf.SurfData);
        if ~isempty(Xline) && ~isempty(Yline) && ~isempty(Zline)
            if ischar(color)
                line(Xline,Yline,Zline,'LineWidth',1.5,'Color',color(s));
            else
                line(Xline,Yline,Zline,'LineWidth',1.5,'Color',color(s,:));
            end
        end
    end
    waitbar(i/ns,h)
end
close(h);
warning('on','MATLAB:conversionToLogical');
warning('on','MATLAB:divideByZero');