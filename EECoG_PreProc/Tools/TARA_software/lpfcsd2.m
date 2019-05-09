function [x, f, cost, lam0, lam1] = lpfcsd2(y, d, fc, theta, sigma, pen, r0, r1, Nit, x_init)
% [x, f, cost] = lpfcsd2(y, d, fc, theta, sigma, pen, r0, r1, Nit)
% Simultaneous low-pass filtering and compound sparsity denoising 
%
% This is a version of 'lpfcsd' that takes (theta, sigma) as input parameters.
% Program 'lpfcsd' takes (lam0, lam1) as input parameters.
%
% INPUT
%   y - raw data
%   d - filter order parameter (d = 1, 2, 3)
%   fc - cut-off frequency (normalized frequency, 0 < fc < 0.5)
%   theta - shape parameter (0 < theta < 1)
%   sigma - noise or psuedo-noise parameter (sigma > 0)
%   pen - penalty function ('L1', 'log', or 'atan')
%   r0, r1 - non-convexity parameters, ignored for 'L1' penalty (0 < r < 1)
%   Nit - number of iterations
%
% OUTPUT
%   x - CSD component
%   f - LPF component
%   cost - cost function history
%
% Use lpfcsd2(..., x_init) to specify initial x.
%
% Use [x, f, cost, lam0, lam1] = lpfcsd2(...) to obtain lambda values.

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
[A, B, B1] = BAfilt(d, fc, N);
Hy = B*(A\y);
BTB = B'*B;
b = BTB*(A\y);
e = ones(N, 1);
D = spdiags([-e, e], [0 1], N-1, N);

trunc = @(x) x(1:end-1);

imp = zeros(N, 1);
imp( round(N/2) ) = 1;

h = B * (A \ imp);
h_norm = sqrt( sum( abs( h ).^2 ) );     % norm of filter

h1 = B1 * trunc(A \ imp);
h1_norm = sqrt( sum( abs( h1 ).^2 ) );     % norm of filter

hh = A' \ (BTB * (A \ imp));
hh_norm = sqrt( sum( abs( hh ).^2 ) );     % norm of filter H^T H

hh1 = B1 * trunc(A \ imp);
hh1 = A' \ (B'* hh1);
hh1_norm = sqrt( sum( abs( hh1 ).^2 ) );     % norm of filter H1^T H

lam0A = 3 * sigma * hh_norm;
lam1A = 3 * sigma * hh1_norm;

lam0 = theta * lam0A;
lam1 = (1-theta) * lam1A;

a0 = r0 * h_norm^2 / lam0A;
a1 = r1 * h1_norm^2 / lam1A;


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

