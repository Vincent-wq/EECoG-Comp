% Generate the complex normal distribution with structured covariance and
% precisioin matrix.
% Reference: Pourahmadi, M. (2013). High-dimensional covariance estimation: 
%            with high-dimensional data (Vol. 882). John Wiley & Sons.
% by Vincent,Pedro @ 2018.12.12

function [PM, COV, Data, nz]=PMgen(p, density, D, n, C, SHOW)
% input:    p: number of variables;
%     density: sparsity of the precision matrix (PM), p*p*density nz entries;
%        Diag: the eigen value of the covariance matrix(COV);
%           n: number of data points generated;
%           C: complex flag;
%        SHOW: show the figure of PM and COV.

% output:  PM: population PM;
%         COV: population COV;
%        Data: p*n data matrix from population with COV and PM.
%          nz: number of non-zero elements.

% fix rand generator 
randn('state',1);
% controling the diag, default identity.
if D~=[]
    DMat = diag(D);
    DMatInv = diag(1./D);
	DMatInvSqrt = diag(1./sqrt(D));
else
    DMat = diag(ones(p,1));
    DMatInv = DMat;
	DMatInvSqrt = DMat;
end
% create sparse triagular matrix with about half spartisty. 
Tr = tril(sprandn(p, p, density));
% for complex and real case
if C   
    Ti = sprandn(Tr);
    T = Tr +1i*Ti;
else
    T=Tr;
end
% change the diagals with ones
T = spdiags(ones(p,1),0,T); 
Tinv = inv(T);
% Obtain the Precision Matrix and Covariance Matrix
PM =T'*DMatInv*T;
COV=Tinv*DMatInv*Tinv';
disp('Sparsity for Cov and PM are: ')
% precise sparsity
nz = [nnz(COV)/(p*p) nnz(PM)/(p*p)]
% display the results
if SHOW
    figure
    subplot(131)
    imagesc(abs(full(COV)))
    colorbar
    title(['COV with sparsity = ' num2str(nz(1))])
    subplot(132)
    imagesc(abs(full(PM)))
    colorbar
    title([' PM with sparsity = ' num2str(nz(2))])
	subplot(133)
    imagesc(abs(full(PM*COV)))
    colorbar
	caxis([0 1])
    title('PM * COV')
end
% Generate Gaussian data from PM
if n~=0
    if C
        Data = Tinv*DMatInvSqrt*(randn(p,n)+1i*randn(p,n));
    else
        Data = Tinv*DMatInvSqrt*randn(p,n);
    end
else
    Data = [];
end
end