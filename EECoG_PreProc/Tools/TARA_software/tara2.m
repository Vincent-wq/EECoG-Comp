function [x1, x2, f, cost, u1, u2, lam0, lam1, lam2, a0, a1, a2] = tara2(y, d, fc, theta, beta, sigma, pen, r0, r1, r2, Nit, u1_init, u2_init)
% [x1, x2, f, cost, u1, u2] = tara2(y, d, fc, theta, beta, sigma, pen, r0, r1, r2, Nit)
% Transient Artifact Reduction Algorithm (TARA)
% Note, this version takes (theta, beta, sigma) as input parameters. 
% The version 'tara' takes (lam0, lam1, lam2) as input parameters.
%
% INPUT
%   y - raw data
%   d - filter order parameter (d = 1, 2, 3)
%   fc - cut-off frequency (normalized frequency, 0 < fc < 0.5)
%   theta, beta - shape parameters (0 < theta < 1, beta > 1)
%   sigma - noise or psuedo-noise parameter (sigma > 0)
%   pen - penalty function ('L1', 'log', or 'atan')
%   r0, r1, r2 - non-convexity parameters, ignored for 'L1' penalty (0 < r < 1)
%   Nit - number of iterations
%
% OUTPUT
%   x1 - sparse signal with sparse derivative
%   x2 - signal with sparse derivative
%   cost - cost function history
%
% Use tara2(..., u1_init, u2_init) to specify initial vectors

y = y(:);
N = length(y);
cost = zeros(1, Nit);

% Use smoothed penalty functions
EPS = 1E-10;
switch pen
    case 'L1'
        phi = @(x, a) sqrt(x.^2 + EPS);
        psi = @(x, a) sqrt(x.^2 + EPS);
        a0 = 0; a1 = 0; a2 = 0;
    case 'log'
        phi = @(x, a) (1/a) * log(1 + a*sqrt(x.^2 + EPS));
        psi = @(x, a) sqrt(x.^2 + EPS) .* (1 + a*sqrt(x.^2 + EPS));
    case 'atan'
        phi = @(x, a) 2./(a*sqrt(3)) .* (atan((2*a.*sqrt(x.^2 + EPS)+1)/sqrt(3)) - pi/6);
        psi = @(x, a) sqrt(x.^2 + EPS) .* (1 + a.*sqrt(x.^2 + EPS) + a.^2.*(x.^2 + EPS));
    otherwise
        disp('Error: penalty must be L1, log, or atan')
        x1 = []; x2 = []; f = []; cost = [];
        return
end  

[A, B, B1] = BAfilt(d,  c, N);
A1 = A(1:N-1, 1:N-1);
e = ones(N, 1);
D = spdiags([-e, e], 0:1, N-1, N);
BTB = B'*B;
B1TB1 = B1'*B1;

% Parameters
trunc = @(x) x(1:end-1);
imp = zeros(N, 1);
imp( round(N/2) ) = 1;


    h = B * (A \ imp);
    h_norm = sqrt( sum( abs( h ).^2 ) );     % norm of filter H
    
    h1 = B1 * trunc(A \ imp);
    h1_norm = sqrt( sum( abs( h1 ).^2 ) );     % norm of filter H1
    
%     hh = A' \ (B' * (B * (A \ imp)));
%     hh_norm = sqrt( sum( abs( hh ).^2 ) );     % norm of filter H^T H

    hh = A' \ (B' * (B * (A \ imp)));
    hh_norm = sqrt( sum( abs( hh ).^2 ) );     % norm of filter H^T H

    hh1 = B1 * trunc(A \ imp);
    hh1 = A' \ (B'* hh1);
    hh1_norm = sqrt( sum( abs( hh1 ).^2 ) );     % norm of filter H1^T H



% hh1 = A' \ (B'* (B1 * trunc(A \ imp)));
% hh1_norm = sqrt( sum( abs( hh1 ).^2 ) );     % norm of filter H1^T H

lam0A = 3 * sigma * hh_norm;
lam1A = 3 * sigma * hh1_norm;

lam0 = lam0A * theta;
lam1 = lam1A * (1 - theta);
lam2 = beta * lam1A;

a0 = r0 * h_norm^2 / lam0;
a1 = r1 * h1_norm^2 / lam1;
a2 = r2 * h1_norm^2 / lam2;


Hy = B * ( A \ y );
y1 = B' * Hy;
y2 = B1' * Hy;

% initialization
if exist('u1_init', 'var')                       
    u1 = u1_init;
else
    u1 = zeros(N, 1);
end
if exist('u2_init', 'var')
    u2 = u2_init;
else
    u2 = zeros(N-1, 1);
end
x1 = A * u1;
Dx2 = A1 * u2;

for i = 1:Nit       
    Lam0 = spdiags( lam0./psi(x1, a0), 0, N, N);
    Lam1 = spdiags( lam1./psi(D*x1, a1), 0, N-1, N-1);
    Lam2 = spdiags( lam2./psi(Dx2, a2), 0, N-1, N-1);
         
    Q1 = 2*BTB + A' * ( Lam0 + D' * Lam1 * D ) * A;
    Q2 = 2*B1TB1 + A1' * Lam2 * A1;

    g = B * u1 - B1 * u2;
    u1 = Q1 \ (y1 + B' * g);
    u2 = Q2 \ (y2 - B1' * g);    
    
    x1 = A * u1;
    Dx2 = A1 * u2;    
    cost(i) = 0.5 * sum((Hy - B*u1 - B1*u2).^2) + ...
        lam0 * sum(phi(x1, a0)) + lam1 * sum(phi(D*x1, a1)) + lam2 * sum(phi(Dx2, a2));
end

x2 = [0; cumsum(Dx2)]; 
bn = nan(d, 1);
f = y - x1 - x2 - [bn; Hy - B*u1 - B1*u2; bn];      % f : low-pass component
