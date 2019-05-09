function [Su] = spd2pd(S)
% correct the very small negative imarinary part due to numerical problem,
% this will make the input S a real semi positive definite matrix.
% 2018.10.29 by Vincent,Pedro
% decompose
[V, D]=eig(S);
% manual correction
d = diag(D);
d(d<0)=0;
% recompose
D=diag(sqrt(d));
Su = V*D*V';
end