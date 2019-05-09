function [S,Data,X] = CSCG(m,q,nblocks,options)
% Pedro Valdes-Sosa, Oct 2017
% Deirel Paz Linares, Oct 2017
% Eduardo Gonzalez-Moreira, Oct 2017
% modified by Vincent for COV output

%% Inputs
% m       - Sample number
% q       - Simulation size 
% nblocks - Number of blocks in the structure
% options - Overlapping Blocks (options = 2) Non-overlapping Blocks (options = 1)
%% Outputs
% S       - Empirical Covariance
% Data    - Samples from random generator
% X       - Partial correlations

if options == 1
    %% Non overlapping blocks
    X            = eye(q); 
    %Infimum and maximum size of blocks
    infsize      = ceil(q/(5*nblocks));
    maxsize      = ceil(q/nblocks) - 1;
    %Cycle by nblocks to fill the matrix X 
    for i = 1:nblocks
        size           = randi([infsize maxsize]);
        blockRe        = 2*(rand(size) - 0.5);
        blockIm        = 2*(rand(size) - 0.5);
        blockRe        = (blockRe + blockRe')/2;
        blockIm        = (blockIm - blockIm')/2;
        block          = blockRe + 1i*blockIm;
        dmin           = min(eig(block));
        if dmin < 0
            block      = block + abs(dmin)*eye(size) + eye(size);
        else
            block      = block - abs(dmin)*eye(size) + eye(size);
        end
        index          = randi([(i-1)*maxsize + 1 i*maxsize],size,1);
        X(index,index) = block; 
    end
end
if options == 2
    %% Overlapping blocks 
    X        = eye(q); 
    %Infimum and maximum size of blocks
    infsize  = ceil(q/(5*nblocks));
    maxsize  = ceil(q/nblocks) - 1;
    %Cycle by nblocks to fill the matrix X 
    for i = 1:nblocks
        size           = randi([infsize maxsize]);
        blockRe        = 2*(rand(size) - 0.5);
        blockIm        = 2*(rand(size) - 0.5);
        blockRe        = (blockRe + blockRe')/2;
        blockIm        = (blockIm - blockIm')/2;
        block          = blockRe + 1i*blockIm;
        if i == nblocks
            index          = randi([(i-1)*maxsize + 1 i*maxsize],size,1);
        else
            index          = randi([(i-1)*maxsize + 1 ceil((i + 1/2)*maxsize)],size,1);
        end
        X(index,index) = block; 
    end
    dmin     = min(eig(X));
    if dmin < 0
        X    = X + abs(dmin)*eye(q) + eye(q);
    else
        X    = X - abs(dmin)*eye(q) + eye(q);
    end    
end
%% Applying isomorphism and generating data
W           = eye(q)/X;
W           = (W + W')/2;
Wisomph     = [real(W) -imag(W); imag(W) real(W)];
Data_isomph = mvnrnd(zeros(1,2*q),Wisomph,m);
DataRe      = Data_isomph(:,1:q);
DataIm      = Data_isomph(:,q+1:2*q);
Data        = DataRe + 1i*DataIm;
S           = (1/m)*(Data'*Data);