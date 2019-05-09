function [H]=calSourLap(S, L, lambda)
% calculate hat matrix from ecog voltage to cortex Laplacian. 
% V 0.00 By Vincent&Pedro
[nChan, ~] = size(S);
L2 = L*L;
S1 = S/L2;
S2 = S1*S1';
H = L\S1'/(S2+lambda*eye(nChan));
end  