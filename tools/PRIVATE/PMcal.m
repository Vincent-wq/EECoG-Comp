function [ThetaJJ] = PMcal(Svv, nS, maxiter11)
% calculate the precision matrix from the cross spectrum matrix Svv;
%       Svv: cross spectrum matrix;
%        nS: number of cortical generators;
% maxiter11: maximum number of iterations.
% By Vincent, Deriel, Pedro @ 2018.11.11

[m,~]= size(Svv);
rho            = sqrt(log(nS)/m);                                          % Regularization parameter
rho_diag       = 0;                                                        % Regularization mask diagonal
rho_ndiag      = 1;                                                        % Regularization mask nondiagonal
lambda2        = (rho*(rho_diag*eye(m)+rho_ndiag*(ones(m)-eye(m)))).^2;    % Regularization mask squared
% maxiter11      = 30;
%% Source Graphical Model Local Quadratic Approximation
[ThetaJJ]      = sggm_lqa(Svv, m, lambda2, rho, maxiter11);                % Source Precision Matrix estimator
% if ifCorrect
%     ThetaJJ    = 2*ThetaJJ - ThetaJJ*Svv*ThetaJJ;                        % Unbiased Source Precision Matrix estimator
% end
%% Source Graphical Model Local Quadratic Approximation            % not standard method
function [PM] = sggm_lqa(EC,m,lambda2,rho,maxiter)
%% sggm-lqa initial values
m2          = m^2;                                                   
m12         = m^(1/2);                                               
nu          = 2.5e0;
EC          = EC + 1e-4*max(abs(diag(EC)))*eye(length(EC));     % update by diag++, like ridge regression
ECinv       = pinv(EC);                                         % sudo inverse
ECinv_ph    = exp(1i*angle(ECinv));                             % 
idx         = (lambda2  > 0);                                   
idx0        = (lambda2 == 0);                                   % diag zeros
gamma       = zeros(length(lambda2));                           
gamma(idx0) = rho*m*abs(ECinv(idx0)).^2;                        % initial. diagnals
PM          = ECinv;
for cont = 1:maxiter
    %% Estimation of variances Gamma of Gaussian Mixtures prior
    DET           = 1 + 4*m2*lambda2(idx).*abs(PM(idx)).^2;
    gamma(idx)    = (sqrt(DET) - 1)./(2*m*lambda2(idx));
    ninf          = max(gamma(:));
    gamman        = gamma/ninf;
    %%
    %% Standarization of the Empirical Covariance Matrix                   % where does this come from?
    st_factor1    = ninf^(-1/2)*ECinv.*gamman.^(1/2)*(nu*m12) + (1/(m12))*ECinv_ph; % Consistent factor1
    st_factor2    = (1 + gamman*(nu*m12));                                          % Consitent factor2
    ECst_inv      = st_factor1./st_factor2;                                         % Inverted Standard Empirical Covariance Matrix
    ECst          = pinv(ECst_inv);                                                 % Standard Empirical Covariance Matrix
    %%
    %% Standard Precision Matrix estimator
    [U,D,V]       = svd(ECst);
    D             = abs(diag(D));
    PMst          = (1/2)*U*diag(sqrt(D.^2 + 4) - D)*V';                 % why 4?
    PMst          = (PMst + PMst')/2;                                    % 
    %%
    %% Unstandarized Precision Matrix estimator
    PM            = gamma.^(1/2).*PMst;
    %%
end
end

end
