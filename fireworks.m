function [x_est] = fireworks(nPoints, x_pred, x_est_GTRS, a_i, d_i, safety_distance, Border)
    % Computing the angle of movement
    theta = atan2(x_pred(2) - x_est_GTRS(2,end), x_pred(1) - x_est_GTRS(1,end)); % Computing the angle of movement

    % Center of the ellipse
    center = x_pred(1:2,end);

    % Minor axis length - half of the distance bewteen x_pred
    % and current x_est
    r_max = norm(x_pred(1:2) - x_est_GTRS(:,end)) / 2;
    
    % Major axis length - half of the distance bewteen x_pred
    % and prior x_est
    d = norm(x_pred(1:2, end) - x_est_GTRS(:,end-1)) / 2;
    % In the second iteration x_est and x_pred are equal
    if d ~= 0
        %------------------------------------------------------------%
        %- Generate points inside ellipse using a circle with radius %
        % 1 entered in [0,0] that is then transformed                %
        % into the desierd ellipse                                  -%
        %------------------------------------------------------------%
        % Generate points inside circle in polar coordinates
        angle = 2 * pi * rand(1,nPoints); % [0, 2 * pi]
        radius = sqrt(rand(1, nPoints)); % sqrt for uniform distribution  

        % Convert generated points from polar coordinates to Cartesian
        x_circle = radius .* cos(angle);
        y_circle = radius .* sin(angle);   

        % Create elipse from circle
        % Scale
        x_circle_scaled = (d/2) * x_circle;
        y_circle_scaled = (r_max/2) * y_circle;

        % Rotate
        x_circle_scaled_rotated =  x_circle_scaled * cos(theta) - y_circle_scaled * sin(theta);
        y_circle_scaled_rotated =  x_circle_scaled * sin(theta) + y_circle_scaled * cos(theta);            
        
        % Center circle in x_pred
        x_circle_scaled_rotated_centered = x_circle_scaled_rotated + center(1);
        y_circle_scaled_rotated_centered = y_circle_scaled_rotated + center(2);
        
        % Combine generated points into a single vector
        points_tot = [x_circle_scaled_rotated_centered; y_circle_scaled_rotated_centered];
        points_tot(:,end + 1) = x_est_GTRS(:,end);

        % Analize points in Elipse and remove the ones outside the
        % scenario (in the case that part of the elipse is outside
        % the scenario, for instace when the estimate is near the borders
        indexes_to_remove = [];
        for index = 1 : 1 : size(points_tot, 2) - 1
            if points_tot(1, index) < 0 || points_tot(1, index) > Border || points_tot(2, index) < 0 || points_tot(2, index) > Border
                indexes_to_remove = [indexes_to_remove, index];
            end
        end
        points_tot(:,indexes_to_remove) = [];

        %--------------------------------------------%
        %- Plot points creation -> Circle to Elipse -%
        %--------------------------------------------%
        % figure;
        % subplot(2,2,1);
        % plot(x_circle, y_circle, 'o');
        % axis equal; grid on;
        % title('1) Circle with radius 1');
        % subplot(2,2,2);
        % plot(x_circle_scaled , y_circle_scaled, 'o');
        % axis equal; grid on;
        % title('2) Circle Scaled');
        % subplot(2,2,3);
        % plot(x_circle_scaled_rotated, y_circle_scaled_rotated, 'o');
        % axis equal; grid on;
        % title('3) Circle Scaled and Rotated');
        % subplot(2,2,4);
        % plot(x_circle_scaled_rotated_centered, y_circle_scaled_rotated_centered, 'o');
        % axis equal; grid on;
        % title('4) Circle Scaled, Rotated an Centered');
       
        %--------------------------%
        %- Plot ellipse for debug -%
        %--------------------------%
        % figure
        % % Create border of the ellipse
        % tt = 0 : pi/100 : 2 * pi;
        % x_ellipse = center(1) + d/2 * cos(tt) * cos(theta) - r_max/2 * sin(tt) * sin(theta);
        % y_ellipse = center(2) + r_max/2 * sin(tt) * cos(theta) + d/2 * cos(tt) * sin(theta);
        % plot(x_ellipse, y_ellipse, 'b', 'LineWidth', 1.5)
        % hold on
        % % Plot d
        % xm = [-d/2 d/2]*cos(theta) - [0 0]*sin(theta) + center(1);
        % ym = [-d/2 d/2]*sin(theta) + [0 0]*cos(theta) + center(2); 
        % plot(xm, ym, 'r', 'LineWidth', 2)
        % % Plot r_max
        % xn = [0 0]*cos(theta) - [-r_max/2 r_max/2]*sin(theta) + center(1);
        % yn = [0 0]*sin(theta) + [-r_max/2 r_max/2]*cos(theta) + center(2);
        % plot(xn, yn, 'g', 'LineWidth', 2)
        % % Plot points
        % scatter(points_tot(1,:), points_tot(2,:), 10, 'r', 'filled')

        % After Fireworks use ML to find the min value of the nPoints
        x_est = MaximumLikelihood(points_tot, a_i, d_i);
    else
        x_est = x_est_GTRS(:,end);
    end
end

function x_est = MaximumLikelihood(points_tot, a_i, d_i)
    % This function implements the Maximum Likelihood function
    N = size(a_i,2); % Number of reference points
    min_ml_value = 10^9; % Minimum value for the ML function
    min_pos = [1000; 1000]; % Position of the minimum value
    % Foreach point inside the elipse compute ML value
    for point = 1 : 1 : size(points_tot,2)
        x = [points_tot(1,point); points_tot(2,point)];
        sum_ml_values = 0;
        % Maximum Likelihood function
        for j = 1 : 1 : N % For each reference point compute ML value
            new_ml_value = ((d_i(j) - norm(x - a_i(:,j)))^2);
            sum_ml_values = sum_ml_values + new_ml_value;
        end
        % Update minimum value and point position
        if sum_ml_values < min_ml_value
            min_ml_value = sum_ml_values;
            min_pos = x;
        end
    end
    % Estimated position
    x_est = min_pos;
end