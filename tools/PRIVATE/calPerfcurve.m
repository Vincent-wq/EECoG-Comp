function [X,Y,T,AUC] = calPerfcurve(compCell,GT,nSample,nBoot)
% calculate the perfcurve for comparison of estimated PMs against PM.
% By Vincent,Pedro @ 2019.1.15

if ~isempty(nBoot); doBoot = 1; else; doBoot =0; end
[dim1, dim2] = size(compCell);
[GTTri, ~] = getTri(GT,1); 
maskGT = (GTTri~=0);

if sum(maskGT)/length(maskGT)>0.6
    [GT] = thresPM(GT, [], 0.05, nSample, 2);
    [GTTri, ~] = getTri(GT,1);
    maskGT = (GTTri~=0);
    disp(['Sparsity: ' num2str(sum(maskGT)/length(maskGT))])
end
    
X = cell(dim1,dim2);   Y = cell(dim1,dim2); 
AUC = cell(dim1,dim2); T = cell(dim1,dim2);

for i1=1:dim1
    for i2 = 1:dim2
        compN = normPM(compCell{i1,i2}, nSample);
        [compNTri, ~] = getTri(compN,1);
        if doBoot
            [X{i1,i2},Y{i1,i2},T{i1,i2},AUC{i1,i2}] = perfcurve(maskGT(:),...
            abs(compNTri(:)), 1, 'NBoot',nBoot, 'BootType','cper','Options',...
            statset('UseParallel',true)); 
        else
            [X{i1,i2},Y{i1,i2},T{i1,i2},AUC{i1,i2}] = perfcurve(maskGT(:),...
            abs(compNTri(:)),1);
        end
    end
end
end