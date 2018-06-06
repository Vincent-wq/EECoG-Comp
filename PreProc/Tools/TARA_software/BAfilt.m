function [A, B, B1, D, a, b, b1] = BAfilt(d, fc, N)
% [A, B, B1, D] = BAfilt(d, fc, N)
%
% Banded matrices for zero-phase high-pass filter.
%   A : [N N]
%   B : [N-2d N]
%   B = B1*D
%   B1 : [N-2d N-1]
%   D : [N-1 N]
% The matrices are 'sparse' data type in MATLAB.
%
% INPUT
%   d  : degree of filter is 2d
%   fc : cut-off frequency (normalized frequency, 0 < fc < 0.5)
%   N  : length of signal
%
% [A, B, B1] = ABfilt() also returns B1 such that B = B1 D
% where D is the first-order different matrix (diff())

b1 = [1 -1];
for i = 1:d-1
    b1 = conv(b1, [-1 2 -1]);
end
b = conv(b1, [-1 1]);

omc = 2*pi*fc;
t = ((1-cos(omc))/(1+cos(omc)))^d;

a = 1;
for i = 1:d
    a = conv(a,[1 2 1]);
end
a = b + t*a;
A = spdiags( a(ones(N,1), :) , -d:d, N, N);     % A: Symmetric banded matrix
A1 = spdiags( a(ones(N-1,1), :) , -d:d, N-1, N-1);
B1 = spdiags(b1(ones(N,1), :) , 0:2*d-1, N-2*d, N-1);       % B1: banded matrix

e = ones(N, 1);
D = spdiags([-e e] , 0:1, N-1, N);
B = B1 * D;


% % Verify that B = B1*D
% x = randn(N,1);
% err = B*x - B1*diff(x);
% if max(abs(err(:))) > 1e-5
%     disp('Error: B not equal to B1*D')
% end

