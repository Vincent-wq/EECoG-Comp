function [DSTF] = mk_imFilters(data, K, lambda, FILTER)
% Create inverse solution DFTF.
% By Vincent, 2019.1.2
[n,~] = size(data);    
switch FILTER
    case 'MNE'
        [mkf]=mkfilt_minnorm(K,lambda);
        DSTF = mkf';
    case 'LCMV'
       [mkf]=mkfilt_lcmv(K,double(corr(data')), lambda); 
       DSTF = mkf';
    case 'eLoreta'
        maxIterELoreta=500; tolEloreta = .0001;  % eLoreta
        tmp = lambda^2;
        [DSTF,~] = eLoretaFilter(K, n*tmp, maxIterELoreta, tolEloreta);
end
end