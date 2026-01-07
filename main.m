%===========================================================
% Elliptical Fireworks
% January 2026
% v1.0.1 - Initial Implementation
% © 2026 COPELABS - Universidade Lusófona CUL
%===========================================================

%===========================================================
% TODO:
% 1 - Ver grossura dos obstáculos no plot?
% 2 - Comparar distâncias previstas com as medidas e ver se ficam dentro de 2
% ou 3 STDs das previstas

% 3 - Testar no cenário do WCM

% X_est_GTRS pode não estar dentro de elipse
%===========================================================

%===========================================================
% Comments:
% Quando qq = 2 x_pred é igual a x_est. Por isso elipse não funciona
%===========================================================

clear variables

%-------------------------%
%- Simulation parameters -%
%-------------------------%
M = 1; % Number of targets
N = 8; % Number of reference points (a_i)
K = 50; % Number of measurement samples
MC = 100; % Monte Carlo runs
Border = 10; % Length of volume of interest
sigma_i = 1; % Noise STD in for distance measurements in meters
moving_step = 0.5; % Step used for moving the UAV
nPoints = 1e3; % Number of points inside the elipse
delta = 0.5; % Object bias
std_obstacle = delta / 10; % Object standard deviation
safety_distance = 0.5; % Safety distance to avoid crashing to the walls
stop_threshold = 1; % Distance to reach the destinations

% Reference Points True Location
% First 4 positions are the corners, the others are middle
a_i = [[0; 0], [Border; Border], [0; Border],[Border;0], [Border/2; 0], [Border/2; Border], [0; Border/2], [Border; Border/2]];

% Obstacles True Location [x1, y1; x2, y2]
obstacles(:,:,1) = [0 2.5; 8 2.5];
%obstacles(:,:,2) = [2 5.5; Border 5.5];
obstacles(:,:,2) = [0 8; 8 8];

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
x_destination = readmatrix('Path.txt')';
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
% Measure running time
tic
while (mc - not_feasible) <= MC
    qq = 1; % Target location counter
    ww = 1;
    x_true = [1; 1];
    x_est = zeros(2,1);
    x_est_GTRS = zeros(2,1);
    x_state = zeros(4,1);
    NLOS_identification = [];

    while ww <= N_dest
        RMSE_goal = [];
        while norm(x_destination(1:2,ww) - x_true(1:2,end)) > stop_threshold
            %----------------------------------------------%
            %- Get measurements from the reference points -%
            %----------------------------------------------%
            x = x_true(:, end);
            % Measurements with influence of objects
            [d_i, d_ik, d_i_clean] = getMeasurments(x_true(:,end),a_i, N, K, sigma_i, obstacles, std_obstacle, delta, safety_distance);
            d_i_all = d_i;
            % Measurements without influence of objects
            %d_i = d_i_clean;

            % Clean NLOS measurments
            if ~isempty(NLOS_identification)
                d_i(NLOS_identification) = d_i(NLOS_identification) - delta_i_hat(NLOS_identification);
            end
            
            %----------------------------%
            %- GTRS position estimation -%
            %----------------------------%
            if qq == 1 % In the first iteration, the first $number_of_anchors_to_use are used to estimate the position
                number_of_anchors_to_use = 4;
                [d_i_aux, index] = sort(d_i, 'ascend');
                d_i = d_i_aux(1:number_of_anchors_to_use);
                a_i = a_i(:,index(1:number_of_anchors_to_use));
            end
            d_weight_ij = []; % Average distance between x and a_i and between x and a_j to form weights
            u_ij = []; % Unit vector between a_i and a_j (@ a_i)
            for ii = 1 : 1 : size(a_i,2) - 1
                for jj = ii + 1 : 1 : size(a_i,2)
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
            % First iteraction uses the GTRS estimation only
            if qq == 1
                eigen_values = eig((A'*A)^(1/2) \ D / (A'*A)^(1/2));
                eig_1 = max(eigen_values);
                min_lim = -1/eig_1; % Lower limit for the considered interval
                max_lim = 1e6; % Upper limit for the considered interval
                tol = 1e-3; % Error tolerance
                N_iter = 30; % Maximum number of iterations for bisection
                lambda = bisection_fun(min_lim, max_lim, tol, N_iter, A, D, b, f); % Calling the bisection function
                y_hat = (A' * A + lambda * D + 1e-6 * eye(3)) \ (A' * b - lambda * f); % Adding regularization term to avoind matrix singularity
                x_est(:, qq) = real(y_hat(1:size(x,1),1)); % y_hat = [x^T, norm(x)^2]^T
                if x_est(1,end) < 0
                   x_est(1,end) = safety_distance;
                elseif x_est(1,end) > Border
                   x_est(1,end) = Border - safety_distance;
                elseif x_est(2,end) < 0
                   x_est(2,end) = safety_distance;
                elseif x_est(2,end) > Border
                   x_est(2,end) = Border - safety_distance;
                end
                x_est_GTRS(:, qq) = x_est(:, qq);
                x_state(:,qq) = [x_est(:,qq); 0; 0]; % Initial target estimation obtained by solving the localization problem
                P = eye(4);
                % Reset a_i to default
                a_i = [[0; 0], [Border; Border], [0; Border],[Border;0], [Border/2; 0], [Border/2; Border], [0; Border/2], [Border; Border/2]];
            else
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
                eigen_values = eig((A_track' * A_track)^(1/2) \ D_track / (A_track' * A_track)^(1/2));
                eig_1 = max(eigen_values);
                min_lim = -1/eig_1; % Lower limit for the considered interval
                lambda_track = bisection_fun(min_lim, max_lim, tol, N_iter, A_track, D_track, b_track, f_track); % I am calling the bisection function
                y_hat_track = (A_track' * A_track + lambda_track * D_track + 1e-6 * eye(size(A_track,2))) \ (A_track' * b_track - lambda_track * f_track); % Adding regularization term to avoind matrix singularity
                % Check if estimation is outside the scenario
                if y_hat_track(1,end) < 0
                    y_hat_track(1,end) = safety_distance;
                elseif y_hat_track(1,end) > Border
                    y_hat_track(1,end) = Border - safety_distance;
                elseif y_hat_track(2,end) < 0
                    y_hat_track(2,end) = safety_distance;
                elseif y_hat_track(2,end) > Border
                    y_hat_track(2,end) = Border - safety_distance;
                end
                x_est_GTRS(:, qq) = y_hat_track(1:size(x,1),1);
                % x_state is the target state (Position and velocity)
                x_state(:,qq) = real(y_hat_track(1:size(x_state,1)));
                P = (x_state(:,qq) - x_state(:,qq-1)) * (x_state(:,qq) - x_state(:,qq-1))';

                %------------------------------------------------%
                %- Position estimate improvement with Fireworks -%
                %------------------------------------------------%
                x_est(:, qq) = fireworks(nPoints,x_pred, x_est_GTRS, a_i, d_i, Border);
                %x_est(:, qq) = x_est_GTRS(:,end);
            end
            
            %--------------------%
            %- Check NLOS links -%
            %--------------------%
            % Estimate delta and std
            delta_i_hat = abs(sum(d_ik - sqrt((x_est(1,end) - a_i(1,:)).^2 + (x_est(2,end) - a_i(2,:)).^2)', 2) / K);
            sigma_hat = sqrt(sum(sum((d_ik - sqrt((x_est(1,end) - a_i(1,:)).^2 + (x_est(2) - a_i(2,:)).^2)' - delta_i_hat).^2, 2) / (N * K - 1)));
            % Predict distances using x_pred and subtract estimated delta
            % and sigma
            d_i_pred = sqrt((x_state(1,end) - a_i(1,:)).^2 + (x_state(2,end) - a_i(2,:)).^2)' - delta_i_hat;
            % Compare measured distances with predicted and check if they
            % are in the range of 2 or 3 sigma_hat
            d_i_compare = abs(d_i_all - d_i_pred);
            %NLOS_identification = find(d_i_compare(:,end) > sigma_hat);
            NLOS_identification = find(d_i_compare > sigma_hat);
            
            % Error between the true and measured distance
            % e_i = abs(sqrt((x_est(1,end) - a_i(1,:)).^2 + (x_est(2,end) - a_i(2,:)).^2 )' - d_i_all);
            % Probability of a link being NLOS
            % p_i = e_i./sum(e_i);
            % Ideal 1/N + sigma
            % NLOS_threshold = 1/N + sigma_i;
            % Comparar dois vetores, se > NLOS, < LOS
            % identification = find(p_i(:,end) > NLOS_threshold);

            %------------------------------%
            %- Move Target using velocity -%
            %------------------------------%
            % Compute velocity using azimute to destination
            azimute = atan2(x_destination(2, ww)- x_est(2,end), x_destination(1, ww) - x_est(1,end));
            uav_velocity = [cos(azimute); sin(azimute)] * moving_step;
            % Compute velocity using the velocity function
            %uav_velocity = velocity(x_est(1:2,qq),x_destination(1:2,ww));
            x_true(1,qq+1) = x_true(1,qq) + uav_velocity(1);
            x_true(2,qq+1) = x_true(2,qq) + uav_velocity(2);

            %-----------------------%
            %- Metrics Calculation -%
            %-----------------------%
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
    plotScenario(obstacles, Border, Border, a_i)
    plot(x_destination(1,:), x_destination(2,:), '-', 'Linewidth', 2.5)
    plot(x_true(1,:), x_true(2,:), '-', 'Linewidth', 2.5)
    mc = mc + 1;
end % while (mc - not_feasible) <= MC

RMSE = [RMSE, sqrt(RMSE_i/MC)];
BIAS = [BIAS, norm( BIAS_i/MC, 1 )];
CDF = [CDF, CDF_i];
% [h,stats] = cdfplot(CDF(:,end)); % Plots the CDF and gives stats
not_feasible_tot = [not_feasible_tot; not_feasible];
toc