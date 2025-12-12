function [d_i, d_i_clean] = getMeasurments(x_true,a_i, N, K, sigma, obstacles, std_obstacle,delta)
    %---------------------------------------------------------------------
    % Compute RSS
    %---------------------------------------------------------------------
    d_ik = zeros(N,K);
    d_ik_clean = zeros(N,K);
    delta_i = 2 * std_obstacle * rand(N,K) + (delta - std_obstacle) * ones(N,K);
    for n = 1 : N
        number_of_interscetions = 0;
        % Compute the number of obsctales to each anchor
        for j = 1 : size(obstacles,3)
            number_of_interscetions = number_of_interscetions + compute_intersections(x_true, a_i(:,n), obstacles(:,:,j));
        end
        d_ik(n,:) = sqrt((x_true(1) - a_i(1,n)).^2 + (x_true(2) - a_i(2,n)).^2)' .* ones(1,K) + (delta_i(n,:) * number_of_interscetions) + sigma * randn(1,K);
        % d_i without the influence of obstacles
        d_ik_clean(n,:) = d_ik(n,:) - (delta_i(n,:) * number_of_interscetions);
    end
    % Median of the distances
    d_i = median(d_ik,2);
    d_i_clean = median(d_ik_clean,2);
end

function has_intersection = compute_intersections(x_true, a_i, obstacles)
    %---------------------------------------------------------------------
    % Obstacle Intersection Simulation
    %---------------------------------------------------------------------
    x1 = a_i(1,1);
    y1 = a_i(2,1);
    x2 = x_true(1);
    y2 = x_true(2);
    x3 = obstacles(1,1);
    y3 = obstacles(1,2);
    x4 = obstacles(2,1);
    y4 = obstacles(2,2);
    t = (((x1-x3)*(y3-y4)) - ((y1-y3)*(x3-x4))) / ((x1-x2)*(y3-y4)-(y1-y2)*(x3-x4));
    u = -(((x1-x2)*(y1-y3)) - ((y1-y2)*(x1-x3))) / ((x1-x2)*(y3-y4)-(y1-y2)*(x3-x4));
    if (t >= 0 && t <=1) && (u >= 0 && u <=1)
        has_intersection = 1;
        % may be usefull in the future
        %intersection_point = [x1 + t*(x2-x1), y1 + t*(y2-y1)];
    else
        has_intersection = 0;
    end
end