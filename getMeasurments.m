function [d_i, d_i_clean] = getMeasurments(x_true, a_i, N, K, sigma, obstacles, std_obstacle, delta, safety_distance)
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
        % Measurments without the influence of obstacles
        d_ik_clean(n,:) = d_ik(n,:) - (delta_i(n,:) * number_of_interscetions);
    end
    % Median of the measured distances
    d_i = median(d_ik,2);
    d_i_clean = median(d_ik_clean,2);
    %---------------------------------------------------------------------
    % Max distance measurment (diagonal between the current anchor and the
    % furthest one
    %---------------------------------------------------------------------
    maxDiagonalCornerAnchors = sqrt((a_i(1,1) - a_i(1,2)).^2 + (a_i(2,1) - a_i(2,2)).^2);
    maxDiagonalMiddleAnchors = sqrt((a_i(1,2) - a_i(1,5)).^2 + (a_i(2,2) - a_i(2,5)).^2);
    % Check if any of the measurments is bigger than the diagonals
    for i = 1 : 1 : 4 % Border Anchors
        if d_i(i) > maxDiagonalCornerAnchors - safety_distance
            d_i(i) = maxDiagonalCornerAnchors - safety_distance;
        end
    end
    for i = 5 : 1 : N % Middle Anchors
        if d_i(i) > maxDiagonalMiddleAnchors - safety_distance
            d_i(i) = maxDiagonalMiddleAnchors- safety_distance;
        end
    end
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
    else
        has_intersection = 0;
    end
end