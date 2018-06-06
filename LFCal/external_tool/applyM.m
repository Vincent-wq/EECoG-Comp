function x = applyM(x,M)

S = size(M);
if all(S == [3 3])
    x = x*M';
elseif all(S == [4 4])
    x = [x ones(size(x,1),1)]*M';
    x(:,4) = [];
else
    error('incorrect dimensions of the transformation matrix');
end