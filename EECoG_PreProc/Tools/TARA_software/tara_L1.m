function [x1, x2, f, cost] = tara_L1(y, d, fc, lam0, lam1, lam2, Nit)
% [x1, x2, f, cost] = tara_L1(y, d, fc, lam0, lam1,lam2, Nit)
% Transient Artifact Reduction Algorithm (TARA)
% with L1 norm penalty
%
% INPUT
%   y - raw data
%   d - filter order parameter (d = 1, 2, 3)
%   fc - cut-off frequency (normalized frequency, 0 < fc < 0.5)
%   lam0, lam1 - regularization parameter for x1 and diff(x1)
%   lam2 - regularization parameter for diff(x2)
%   Nit - number of iterations
%
% OUTPUT
%   x1 - sparse signal with sparse derivative
%   x2 - signal with sparse derivative
%   cost - cost function history

y = y(:);
N = length(y);
cost = zeros(1, Nit);
EPS = 1E-10;
phi = @(x) sqrt(x.^2 + EPS);

[A, B, B1] = BAfilt(d, fc, N);
A1 = A(1:N-1, 1:N-1);
e = ones(N, 1);
D = spdiags([-e, e], 0:1, N-1, N);
BTB = B'*B;
B1TB1 = B1'*B1;

Hy = B * ( A \ y );
y1 = B' * Hy;
y2 = B1' * Hy;

u1 = zeros(N, 1);
u2 = zeros(N-1, 1);
x1 = zeros(N, 1);
Dx2 = zeros(N-1, 1);

for i = 1:Nit       
    Lam0 = spdiags( lam0./phi(x1), 0, N, N);
    Lam1 = spdiags( lam1./phi(D*x1), 0, N-1, N-1);
    Lam2 = spdiags( lam2./phi(Dx2), 0, N-1, N-1);
         
    Q1 = 2*BTB + A' * ( Lam0 + D' * Lam1 * D ) * A;
    Q2 = 2*B1TB1 + A1' * Lam2 * A1;

    g = B * u1 - B1 * u2;
    u1 = Q1 \ (y1 + B' * g);
    u2 = Q2 \ (y2 - B1' * g);    
    
    x1 = A * u1;
    Dx2 = A1 * u2;    
    cost(i) = 0.5 * sum((Hy - B*u1 - B1*u2).^2) + ...
        lam0 * sum(phi(x1)) + lam1 * sum(phi(D*x1)) + lam2 * sum(phi(Dx2));
end

x2 = [0; cumsum(Dx2)]; 
bn = nan(d, 1);
f = y - x1 - x2 - [bn; Hy - B*u1 - B1*u2; bn];      % f : low-pass component
