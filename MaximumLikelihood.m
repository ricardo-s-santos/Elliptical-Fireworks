function x_est = MaximumLikelihood(points_tot, sigma_i, a_i, d_i)
    % This function implements the Maximum Likelihood
    N = size(a_i,2);
    min_ml_value = 10^9;
    min_pos = [1000; 1000];
    for particle = 1 : 1 : size(points_tot,2)
        % Eq. (2) Maximum Likelihood
        x = [points_tot(1,particle); points_tot(2,particle)];
        sum_ml_values = 0;
        for j = 1 : 1 : N
            new_ml_value = (1 / (2 * (sigma_i^2))) * ((d_i(j) - norm(x - a_i(:,j)))^2);
            sum_ml_values = sum_ml_values + new_ml_value;
        end
        if sum_ml_values < min_ml_value
            min_ml_value = sum_ml_values;
            min_pos = x;
        end
    end
    x_est = min_pos;
end