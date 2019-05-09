function [DSTF] = cross_nonovgrouped_elasso_ssbl(Svv,LeadFields,Nsamp,nonovgroups)
% Elitist Lasso_Sparse Bayesian Learning
%% Loading Dimensions
Nv         = size(LeadFields{1},2);                                     % Number of nodes
Ngroups    = length(nonovgroups);                                       % Number of spatial groups
Nt         = size(Svv{1},3);                                            % Number of Time points or Frequency bins
Nsubj      = size(Svv,2);                                               % Number of Subjects
sigma      = ones(Nv,1);                                                % Prior Variances
sigma_post = ones(Nv,Nsubj);                                            % Posterior Variances
miu_cal    = ones(Nv,Nt,Nsubj);                                         % Inverse Solution Calibration
miu        = ones(Nv,Nt,Nsubj);                                         % Inverse Solution
miu2group  = ones(Ngroups,1);                                           % Squered Inverse Solution by groups
DSTF       = cell(1,Nsubj);                                             % Data to Source Transfer Function
h          = ones(Nv,1);                                                % Hypeparameter for Gamma
delta      = ones(Nv,1);                                                % Hypeparameter of spatial correlations
%% User defined parameters
maxiter1   = 15;                                                        % Number of Iterations of outer cycle
maxiter11  = 15;                                                        % Number of Iterations of inner cycle
ealpha     = 1E-8;                                                      % Rate parameter of the gamma pdf of alpha (higher -> smoother)
delta_0    = 1E-2;                                                      % Minimum value of delta's
beta       = 1;                                                         % Data variance
w          = ones(Ngroups,1);                                           % Vector of parameters' weights
W          = repmat(transpose(w),Ngroups,1)-diag(w);                    % Parameters to Deltas' Transformation matrix
%% Initialization of parameters
for cont4 = 1:Nsubj
    K_tmp            = LeadFields{cont4};
    Kt_tmp           = transpose(K_tmp);
    Ne_tmp           = size(K_tmp,1);
    Ine_tmp          = spdiags(ones(Ne_tmp,1),0,Ne_tmp,Ne_tmp);         % Identity matrix with Ne size
    scale_K_tmp      = sqrt(trace(K_tmp*Kt_tmp)/Ne_tmp);                % Lead Field Scale Norm inf
    K_tmp            = K_tmp/scale_K_tmp;                               % Lead Field scaling
    Kt_tmp           = Kt_tmp/scale_K_tmp;
    %% Calibration Inverse Solution
    sigmaKt          = spdiags(sigma,0,Nv,Nv)*Kt_tmp;
    sigma_post1      = sigmaKt/(K_tmp*sigmaKt+beta*Ine_tmp);
    sigma_post2      = K_tmp*spdiags(sigma,0,Nv,Nv);
    for jj=1:Nv
        sigma_post(jj,cont4) = sigma(jj)-sigma_post1(jj,:)*sigma_post2(:,jj);
    end
    %% Compute 'miu_cal' for all slices of 'V'
    DSTF_subj        = (1/beta).*(sigmaKt - sigma_post1*(sigma_post2*Kt_tmp));
    DSTF{cont4}      = DSTF_subj;
    for cont3 = 1:Nt
        miu_cal(:,cont3,cont4) = diag(DSTF_subj*squeeze(Svv{cont4}(:,:,cont3))*DSTF_subj');
    end
    Ne{cont4}        = Ne_tmp;
    Ine{cont4}       = Ine_tmp;
    K{cont4}         = K_tmp;
    Kt{cont4}        = Kt_tmp;
    KtK{cont4}       = Kt_tmp*K_tmp;
    scale_K{cont4}   = scale_K_tmp;
    scale_V{cont4}   = (sum(mean(reshape(abs(miu_cal(:,:,cont4)),Nv,Nt),2))/Nv)/max(sigma_post(:,cont4));
end
%% Data Scaling
for cont4 = 1:Nsubj
    Svv{cont4} = Svv{cont4}/scale_V{cont4};
end
%% Initialization of Hyperparameters
alpha        = 1E-4;                                                                         % Hyperparameter of Elitist LASSO norm 
%% Main Cycle
for cont1 = 1:maxiter1
    for cont11 = 1:maxiter11
        for cont4 = 1:Nsubj
            %% Update Posterior Mean and Covariance matrix
            sigmaKt     = spdiags(sigma,0,Nv,Nv)*Kt{cont4};
            sigma_post1 = sigmaKt/(K{cont4}*sigmaKt+beta*Ine{cont4});
            sigma_post2 = K{cont4}*spdiags(sigma,0,Nv,Nv);
            % Only save the diagonals of the Posterior Covariance
            for jj=1:Nv
                sigma_post(jj,cont4) = sigma(jj)-sigma_post1(jj,:)*sigma_post2(:,jj);
            end
            % Iterative Transference Operator
            DSTF_subj        = (1/beta).*(sigmaKt-sigma_post1*(sigma_post2*Kt{cont4}));
            % Compute 'miu' for all slices of 'V'
            for cont3 = 1:Nt
                miu(:,cont3,cont4) = diag(DSTF_subj*squeeze(Svv{cont4}(:,:,cont3))*DSTF_subj');
            end
        end
        miu2            = sum(reshape(abs(miu),Nv,Nt*Nsubj),2)/(Nt*Nsubj);
        sigma_post_mean = sum(sigma_post,2)/Nsubj;
        %% Update Gammas
        for group = 1:Ngroups
            idx_group         = nonovgroups(group);
            idx_group         = idx_group{1};
            miu2group(group)  = sum(miu2(idx_group));
        end
        delta_group           = W*(Nt*Nsamp*Nsubj*miu2group).^(1/2)+delta_0.*ones(Ngroups,1);
        for group = 1:Ngroups
            idx_group         = nonovgroups(group);
            idx_group         = idx_group{1};
            delta(idx_group)  = delta_group(idx_group);
            h(idx_group)      = sqrt((1./4).^2+alpha^2.*Nt.*Nsamp.*Nsubj.*(w(group).^2).*sum(miu2(idx_group) + sigma_post_mean(idx_group)).*delta_group(idx_group).^2) - 1./4;
        end
        index_h               = find((miu2+sigma)<0);
        h(index_h)            = 0;
        gamma                 = alpha.*delta.^2+h;
        sigma                 = 1/(2*alpha)-delta.^2./(2*gamma);
        sigma_bar             = 2*alpha*sigma;
    end
    %% Update alpha
    index_alpha               = find(sigma_bar>0.005*max(max(sigma_bar)));
    mixed_norm                = (sum((Nt*Nsamp*Nsubj*miu2).^(1/2)))^2;
    c1                        = (Nt*Nsamp*Nsubj*sum((miu2(index_alpha)+sigma_post_mean(index_alpha))./(sigma_bar(index_alpha)))+sum(sum(delta.^2./(1-sigma_bar)))+mixed_norm)/(Nv);
    f_aux                     = @(alpha_aux) c1+ealpha*Nt*Nsamp*Nsubj-Nt*Nsamp*Nsubj/alpha_aux-(1/2)*(1+length(index_alpha)/(Nv))*(1/alpha_aux)-(1/(Nv))*sum(sum((trascendent_term(alpha_aux*delta(:).^2))'.*(delta(:)).^2));
    alpha                     = fzero(f_aux,[10^(-8) 100000]);
end
for cont4 = 1:Nsubj
DSTF{cont4} = (DSTF{cont4}*scale_V{cont4}/scale_K{cont4});
end
end