%**************************************************************************
% Find the best lambda (penalized parameter) by calculating GCV, optimized 
% for the ECoG projections problems.GCV calculated according to Wahba (1979).
% Problem definition: Vecog = S*L*LecogS+error;
%                     s.t. min||Vecog-S*L*LecogS||^2+lambda||L*L*LecogS||^2
%            GCV = n * (H*y)^2 / trace(H)^2
%              H = I-X\(X'X+n*lambda*I)*X' (same case for ridge regression)
% V 0.00 by Vincent&Pedro @ 2018.11.7
% Version 2018.12.15 updated by Vincent @ Cubda.
%**************************************************************************
function [lambda, LAMBDA, GCV] = LambdaGCV(S, L, Vecog, N,logStart, LogStop, NAME)
%       input:      S: selection matrix;
%                   L: Laplacian on cortex;
%               Vecog: data matrix, y in general model;
% N,logStart, LogStop: define the scope of lambda, LAMBDA = 
%                      logspace(logStart, LogStop, N);
%                NAME: name of result figure to save. 
%       output:   GCV: GCV vector;
%              LAMBDA: lambda vector;
%              lambda: best lambda;
LAMBDA = logspace(logStart, LogStop, N);
[nChan,~] = size(Vecog);
GCV=zeros(1,N);
L2 = L*L; 
S1 = S/L2;
S2 = S1*S1';                                   
for i = 1:N
    H = eye(nChan)-S2/(S2+LAMBDA(i)*eye(nChan));
    GCV(i) = nChan*norm(H*Vecog,'fro').^2/trace(H).^2;
end

[~, indLambda] = min(GCV)
lambda=LAMBDA(indLambda)

figure
loglog(LAMBDA, GCV)
xlabel('\lambda')
ylabel('GCV')
title(['GCV V.S. \lambda, min(\lambda)=',num2str(lambda)])
if ischar(NAME)
    fig=gcf;
    FILEPATH = strcat(NAME,['-lambdaGCVind=',num2str(indLambda)], '.jpg');
    saveas(fig,FILEPATH);
end
end