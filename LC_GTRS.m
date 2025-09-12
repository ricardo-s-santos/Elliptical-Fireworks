clear variables
tic

%-------------------%
%- Simulation part -%
%-------------------%

M = 1; % Number of targets
N = 6; % Number of a_i
K = 10; % Number of measurement samples
Border = 100; % Length of volume of interest
sigma_i = 10; % Noise STD in for distance measurements in meters

MC = 1e4; % Monte Carlo runs
RMSE = []; % Root mean square error
BIAS = []; % Bias of the estimator
CDF = []; % Cumulative distributed function
not_feasible_tot = []; % Total not feasible cases

RMSE_i = 0; % Error in each MC run
BIAS_i = 0; % Bias in each MC run
CDF_i = []; % CDF in each MC run
not_feasible = 0; % In the case the method doesn't work for all cases

counter = 1;
while (counter - not_feasible) <= MC
    
    scenario = 0;
    while scenario == 0
        stoping = 0;
        a_i = rand(2,N) * Border;
%         a_i = [[0; 0], [Border; 0], [Border; Border], [0; Border]];
        x = rand(2,M) * Border;
        %========================
        % ||x_i - x_j|| <= 1   
        %========================
        for i = 1 : N
            for j = 1 : M
                anchor_a = a_i(:,i);
                target_b = x(:,j);
                if norm(anchor_a - target_b) < 1
                    stoping = 1;
                    break
                end
            end
        end
        for i = 1 : (N-1)
            for j = (i+1) : N
                anchor_a = a_i(:,i);
                anchor_b = a_i(:,j);
                if norm(anchor_a - anchor_b) < 1
                    stoping = 1;
                    break
                end
            end
        end
        if stoping == 0
            scenario = 1;
        end
    end
    d_ik = sqrt( (x(1,1) - a_i(1,:)).^2 + (x(2,1) - a_i(2,:)).^2 )' + sigma_i * randn(N,K);

    %-------------------%
    %- Estimation part -%
    %-------------------%
    
    d_i = median(d_ik,2);
    d_weight_ij = []; % Average distance between x and a_i and between x and a_j to form weights
    u_ij = []; % Unit vector between a_i and a_j (@ a_i)
    for ii = 1 : 1 : N-1
        for jj = ii+1 : 1 : N
            u_ij = [u_ij, (a_i(:,jj) - a_i(:,ii))/norm((a_i(:,jj) - a_i(:,ii)))];
            d_weight_ij = [d_weight_ij; ii, jj, ( d_i(ii) + d_i(jj) )/2];
        end
    end
    w_ij = (1./d_weight_ij(:,3))./(sum(1./d_weight_ij(:,3))); % Weights
    A1 = [];
    b1 = [];
    for tt = 1 : 1 : size(u_ij,2)
        ii = d_weight_ij(tt,1);
        jj = d_weight_ij(tt,2);
        A1 = [A1; 2 *(norm(a_i(:,ii) - a_i(:,jj)) * u_ij(:,tt) + a_i(:,ii))', -1];
%         A1 = [A1; 2 *(norm(a_i(:,ii) - a_i(:,jj)) * u_ij(:,tt) - a_i(:,jj))', 1];
%         A1 = [A1; 2 *(a_i(:,ii) + a_i(:,jj))', -2];
        b1 = [b1; norm(a_i(:,ii) - a_i(:,jj))^2 - d_i(jj)^2 + 2 * norm(a_i(:,ii) - a_i(:,jj)) * u_ij(:,tt)' * a_i(:,ii) + norm(a_i(:,ii))^2];
%         b1 = [b1; norm(a_i(:,ii) - a_i(:,jj))^2 + d_i(ii)^2 + 2 * norm(a_i(:,ii) - a_i(:,jj)) * u_ij(:,tt)' * a_i(:,ii) - norm(a_i(:,jj))^2];
%         b1 = [b1; norm(a_i(:,ii) - a_i(:,jj))^2 - d_i(ii)^2 - d_i(jj)^2 + 2 * a_i(:,ii)' * a_i(:,jj)];
    end
    W = diag(sqrt(w_ij')); % Weight matrix
%     W = eye(size(d_weight_ij,1)); % If no weights are employed
    A = W * A1;
    b = W * b1;
    D = eye(size(x,1)+1); D(size(x,1)+1,size(x,1)+1) = 0;
    f = [zeros(size(x,1),1); -1/2];
    eigen_values = eig( (A'*A)^(1/2) \ D / (A'*A)^(1/2) );
    eig_1 = max(eigen_values);
    min_lim = -1/eig_1; % Lower limit for the considered interval
    max_lim = 1e6; % Upper limit for the considered interval
    tol = 1e-3; % Error tolerance
    N_iter = 30; % Maximum number of iterations for bisection
    lambda = bisection_fun(min_lim, max_lim, tol, N_iter, A, D, b, f); % Calling the bisection function
    y_hat = (A' * A + lambda * D + 1e-6 * eye(3)) \ (A' * b - lambda * f); % Adding regularization term to avoind matrix singularity
    x_est = y_hat(1:size(x,1),1); % y_hat = [x^T, norm(x)^2]^T

    %---------------------%
    %- Error Calculation -%
    %---------------------%

    if sum(sum(isnan(x_est))) >= 1 || sum(sum(isinf(x_est))) >= 1
        not_feasible = not_feasible + 1;
    else
        RMSE_i = RMSE_i + (x - x_est)' * (x - x_est);
        BIAS_i = BIAS_i + (x_est - x);
        CDF_i = [CDF_i; norm(x - x_est)/M];
    end
    counter = counter + 1;
end % while (counter - not_feasible) <= MC

RMSE = [RMSE, sqrt(RMSE_i/MC)]
BIAS = [BIAS, norm( BIAS_i/MC, 1 )];
CDF = [CDF, CDF_i];
% [h,stats] = cdfplot(CDF(:,end)); % Plots the CDF and gives stats
not_feasible_tot = [not_feasible_tot; not_feasible];

% fName = strcat('LC_GTRS_2D_var_N_sigma',int2str(sigma_i),'_K',int2str(K),'_B',int2str(Border),'_MC',int2str(MC));
% save(fName)

toc