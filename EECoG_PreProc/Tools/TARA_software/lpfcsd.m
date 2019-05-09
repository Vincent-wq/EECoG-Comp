function [x, f, cost] = lpfcsd(y, d, fc, lam0, lam1, pen, a0, a1, Nit, x_init)
% [x, f, cost] = lpfcsd(y, d, fc, lam0, lam1, pen, a0, a1, Nit)
% Simultaneous low-pass filtering and compound sparsity denoising 
%
% INPUT
%   y - raw data
%   d - filter order parameter (d = 1, 2, 3)
%   fc - cut-off frequency (normalized frequency, 0 < fc < 0.5)
%   lam0, lam1 - regularization parameters for x and diff(x)
%   pen - penalty function ('L1', 'log', or 'atan')
%   a0, a1 - non-convexity parameters (ignored for 'L1' penalty)
%   Nit - number of iterations
%
% OUTPUT
%   x - CSD component
%   f - LPF component
%   cost - cost function history
%
% Use lpfcsd(..., x_init) to specify initial x

% Ivan Selesnick
% NYU Polytechnic School of Engineering, New York, USA
% September 2013
% revised February 2014

% Use smoothed penalty functions
EPS = 1E-10;
switch pen
    case 'L1'
        phi = @(x, a) sqrt(x.^2 + EPS);
        psi = @(x, a) sqrt(x.^2 + EPS);
        a0 = 0;
        a1 = 0;
    case 'log'
        phi = @(x, a) (1/a) * log(1 + a*sqrt(x.^2 + EPS));
        psi = @(x, a) sqrt(x.^2 + EPS) .* (1 + a*sqrt(x.^2 + EPS));
    case 'atan'
        phi = @(x, a) 2./(a*sqrt(3)) .* (atan((2*a.*sqrt(x.^2 + EPS)+1)/sqrt(3)) - pi/6);
        psi = @(x, a) sqrt(x.^2 + EPS) .* (1 + a.*sqrt(x.^2 + EPS) + a.^2.*(x.^2 + EPS));
    otherwise
        disp('Error: penalty must be L1, log, or atan')
        x = []; f = []; cost = [];
        return
end

cost = zeros(1, Nit);       % cost function history
y = y(:);                   % convert to column vector
N = length(y);
[A, B] = BAfilt(d, fc, N);
Hy = B*(A\y);
BTB = B'*B;
b = BTB*(A\y);
e = ones(N, 1);
D = spdiags([-e, e], [0 1], N-1, N);

if exist('x_init', 'var')                       % initialization
    x = x_init;
else
    x = y;
end

for k = 1:Nit
    Lam0 = spdiags( lam0./psi(x, a0), 0, N, N);
    Lam1 = spdiags( lam1./psi(D*x, a1), 0, N-1, N-1);
    Q = BTB + A' * (Lam0 + D'*Lam1*D) * A;
    u = Q \ b;
    x = A * u;      % Note: H * x = B * (A\x) = B * u    
    cost(k) = lam0 * sum(phi(x, a0)) + lam1 * sum(phi(D*x, a1)) + 0.5 * sum(abs(B*u-Hy).^2);
end

bn = nan(d, 1);                     % bn : nan's to extend f to length N
f = y - x - [bn; Hy - B*u; bn];     % f : low-pass component

