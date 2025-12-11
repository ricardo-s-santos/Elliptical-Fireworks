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