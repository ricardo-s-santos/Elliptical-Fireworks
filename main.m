%===========================================================
% Elliptical Fireworks
% September 2025
% v1.0 - Initial Implementation
% © 2025 COPELABS - Universidade Lusófona CUL
%===========================================================
clear variables
tic

%-------------------------%
%- Simulation parameters -%
%-------------------------%
M = 1; % Number of targets
N = 6; % Number of a_i
K = 10; % Number of measurement samples
MC = 100; % Monte Carlo runs
Border = 10; % Length of volume of interest
sigma_i = 0.1; % Noise STD in for distance measurements in meters
moving_step = 0.1; % Step used for moving the UAV
nPoints = 1e3; % Number of points inside the elipse

% Reference points true location
a_i = [[0; 0], [0; Border], [Border/2; 0], [Border/2; Border], [Border; Border],[Border;0]];

% State Matrix
delta_t = 1;
S = eye(4);
S(1, 3) = delta_t;
S(2, 4) = delta_t;

% Covariance Matrix
q = 0.05; % State Process Noise (m^2/s^3)
Q = q * [[(delta_t^3)/3, 0, (delta_t^2)/2, 0];
    [0, (delta_t^3)/3, 0, (delta_t^2)/2];
    [(delta_t^2)/2, 0, delta_t, 0];
    [0, (delta_t^2)/2, 0, delta_t]
    ];

% Read waypoints from file
x_destination = csvread('Path.txt')';
N_dest = size(x_destination,2); % Number of destinations

% Metrics variables
RMSE = []; % Root mean square error
BIAS = []; % Bias of the estimator
CDF = []; % Cumulative distributed function
not_feasible_tot = []; % Total not feasible cases

RMSE_i = 0; % Error in each MC run
BIAS_i = 0; % Bias in each MC run
CDF_i = []; % CDF in each MC run
not_feasible = 0; % In the case the method doesn't work for all cases

mc = 1;
while (mc - not_feasible) <= MC
    qq = 1; % Target location counter
    ww = 1;
    x_true = [1; 1];
    x_est = zeros(2,1);
    x_state = zeros(4,1);

    while ww <= N_dest
        RMSE_goal = [];
        while norm(x_destination(1:2,ww) - x_est(1:2,end)) > 0.1
            %--------------------%
            %- Get measurements -%
            %--------------------%
            x = x_true(:, end);
            d_ik = sqrt((x(1,1) - a_i(1,:)).^2 + (x(2,1) - a_i(2,:)).^2 )' + sigma_i * randn(N,K);
            d_i = median(d_ik,2);
            %-------------------%
            %- Estimation part -%
            %-------------------%
            d_weight_ij = []; % Average distance between x and a_i and between x and a_j to form weights
            u_ij = []; % Unit vector between a_i and a_j (@ a_i)
            for ii = 1 : 1 : N-1
                for jj = ii+1 : 1 : N
                    u_ij = [u_ij, (a_i(:,jj) - a_i(:,ii))/norm((a_i(:,jj) - a_i(:,ii)))];
                    d_weight_ij = [d_weight_ij; ii, jj, (d_i(ii) + d_i(jj))/2];
                end
            end
            w_ij = (1./d_weight_ij(:,3))./(sum(1./d_weight_ij(:,3))); % Weights
            A1 = [];
            b1 = [];
            for tt = 1 : 1 : size(u_ij,2)
                ii = d_weight_ij(tt,1);
                jj = d_weight_ij(tt,2);
                A1 = [A1; 2 *(norm(a_i(:,ii) - a_i(:,jj)) * u_ij(:,tt) + a_i(:,ii))', -1];
                b1 = [b1; norm(a_i(:,ii) - a_i(:,jj))^2 - d_i(jj)^2 + 2 * norm(a_i(:,ii) - a_i(:,jj)) * u_ij(:,tt)' * a_i(:,ii) + norm(a_i(:,ii))^2];
            end
            W = diag(sqrt(w_ij')); % Weight matrix
            A = W * A1;
            b = W * b1;
            D = eye(size(x,1)+1); D(size(x,1)+1,size(x,1)+1) = 0;
            f = [zeros(size(x,1),1); -1/2];
            % Using prediction
            if qq ~= 1
                P_pred = S * P * S' + Q;
                x_pred = S * x_state(:, end);
                A1_track = [];
                b1_track = [];
                for tt = 1 : 1 : size(u_ij,2)
                    ii = d_weight_ij(tt,1);
                    jj = d_weight_ij(tt,2);
                    % Updated A1_track
                    A1_track = [A1_track; 2 *(norm(a_i(:,ii) - a_i(:,jj)) * u_ij(:,tt) + a_i(:,ii))' zeros(1,size(x,1)) -1];
                    b1_track = [b1_track; norm(a_i(:,ii) - a_i(:,jj))^2 - d_i(jj)^2 + 2 * norm(a_i(:,ii) - a_i(:,jj)) * u_ij(:,tt)' * a_i(:,ii) + norm(a_i(:,ii))^2];
                end
                A1_track = [A1_track; P_pred^(-1/2), zeros(size(x_state,1),1)];
                b1_track = [b1_track; P_pred^(-1/2) * x_pred];
                W_track = diag([sqrt(1/2) * w_ij', sqrt(1/8) * ones(1,4)]);
                D_track = zeros(5); D_track(1,1) = 1; D_track(2,2) = 1;
                f_track = [zeros(4,1); -1/2];
                A_track = W_track * A1_track;
                b_track = W_track * b1_track;
            end
            % First iteraction uses the estimation only
            if qq == 1
                eigen_values = eig((A'*A)^(1/2) \ D / (A'*A)^(1/2));
                eig_1 = max(eigen_values);
                min_lim = -1/eig_1; % Lower limit for the considered interval
                max_lim = 1e6; % Upper limit for the considered interval
                tol = 1e-3; % Error tolerance
                N_iter = 30; % Maximum number of iterations for bisection
                lambda = bisection_fun(min_lim, max_lim, tol, N_iter, A, D, b, f); % Calling the bisection function
                y_hat = (A' * A + lambda * D + 1e-6 * eye(3)) \ (A' * b - lambda * f); % Adding regularization term to avoind matrix singularity
                x_est(:, qq) = y_hat(1:size(x,1),1); % y_hat = [x^T, norm(x)^2]^T
                x_state(:,qq) = [x_est(:,qq); moving_step; 0]; % Initial target estimation obtained by solving the localization problem
                P = eye(4);
            else % Use prediction
                eigen_values = eig((A_track'*A_track)^(1/2) \ D_track / (A_track'*A_track)^(1/2));
                eig_1 = max(eigen_values);
                min_lim = -1/eig_1; % Lower limit for the considered interval
                lambda_track = bisection_fun(min_lim, max_lim, tol, N_iter, A_track, D_track, b_track, f_track); % I am calling the bisection function
                y_hat_track = (A_track' * A_track + lambda_track * D_track + 1e-6 * eye(size(A_track,2))) \ (A_track' * b_track - lambda_track * f_track); % Adding regularization term to avoind matrix singularity
                x_est(:, qq) = y_hat_track(1:size(x,1),1);
                x_state(:,qq) = real(y_hat_track(1:size(x_state,1)));
                P = (x_state(:,qq) - x_state(:,qq-1)) * ( x_state(:,qq) - x_state(:,qq-1))';
                %--------------------%
                %-  Fireworks part  -%
                %--------------------%
                % figure
                % hold on
                theta = atan2(x_pred(2) - x_est(2,end), x_pred(1) - x_est(1,end)); % Computing the angle of movement
                center = x_est(:,end); % Center of the ellipse
                r_max = norm(x_est(:,end) - x_est(:,end-1)); % Minor axis length
                d = norm(x_est(:,end) - x_pred(1:2)); % Major axis length
                tt = 0 : pi/100 : 2 * pi;
                x_ellipse = center(1) + d/2 * cos(tt) * cos(theta) - r_max/2 * sin(tt) * sin(theta);
                y_ellipse = center(2) + r_max/2 * sin(tt) * cos(theta) + d/2 * cos(tt) * sin(theta);
                % plot(x_ellipse, y_ellipse, 'b')
                % plot(x_est(1,end), x_est(2,end), 'rx', 'MarkerSize',10)
                % plot(x_est(1,end-1), x_est(2,end-1), 'kx','MarkerSize',10)
                % plot(x_pred(1), x_pred(2), 'ro', 'MarkerSize',10)

                pointsInEllipse = 0;
                % while (pointsInEllipse < nPoints)
                %     points_tot = center + d * rand(2, nPoints) - d/2 * ones(2, nPoints); % Generating random points within a square region
                %     plot(points_tot(1,:), points_tot(2,:),'ys')
                %     inEllipse = ((points_tot(1,:) - center(1)) * cos(theta) + (points_tot(2,:) - center(2)) ...
                %       * sin(theta)).^2/(d/2)^2 + (-(points_tot(1,:)-center(1)) * sin(theta) + (points_tot(2,:) - center(2)) ...
                %       * cos(theta)).^2/(r_max/2)^2 <= 1; % Determining points inside the ellipse
                %     pointsInEllipse = sum(inEllipse);
                %     gg = 0;
                % end
                % plot(points_tot(1,inEllipse), points_tot(2,inEllipse),'g*')
                % After Fireworks use ML to find the min value of the
                % particles

            end
            %------------------------------%
            %- Move Target using velocity -%
            %------------------------------%
            % Compute velocity using azimute to destination
            azimute = atan2(x_destination(2, ww)- x_est(2,qq), x_destination(1, ww) - x_est(1,qq));
            %uav_velocity = [cos(azimute); sin(azimute)] * moving_step;
            % Compute velocity using the velocity function
            uav_velocity = velocity(x_est(1:2,qq),x_destination(1:2,ww));
            x_true(1,qq+1) = x_true(1,qq) + uav_velocity(1);
            x_true(2,qq+1) = x_true(2,qq) + uav_velocity(2);
            %---------------------%
            %- Error Calculation -%
            %---------------------%
            if sum(sum(isnan(x_est))) >= 1 || sum(sum(isinf(x_est))) >= 1
                not_feasible = not_feasible + 1;
            else
                RMSE_i = RMSE_i + (x - x_est(:,qq))' * (x - x_est(:,qq));
                BIAS_i = BIAS_i + (x_est(:,qq) - x);
                CDF_i = [CDF_i; norm(x - x_est(:,qq))/M];
            end
            qq = qq + 1;
        end
        ww = ww + 1;
    end
    figure
    plotScenario(Border,a_i)
    plot(x_destination(1,:), x_destination(2,:), '-', 'Linewidth', 2.5)
    plot(x_est(1,:), x_est(2,:), '-', 'Linewidth', 2.5)
    mc = mc + 1;
end % while (mc - not_feasible) <= MC

RMSE = [RMSE, sqrt(RMSE_i/MC)]
BIAS = [BIAS, norm( BIAS_i/MC, 1 )];
CDF = [CDF, CDF_i];
% [h,stats] = cdfplot(CDF(:,end)); % Plots the CDF and gives stats
not_feasible_tot = [not_feasible_tot; not_feasible];
toc