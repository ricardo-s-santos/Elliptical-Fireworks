function lambda = bisection_fun(min_lim, max_lim, tol, N_iter, A, D, b, f)
    lambda = (min_lim + max_lim)/2;
    fun_val = 10^9;
    num_iter = 1;
    while num_iter <= N_iter && abs(fun_val) > tol && abs(min_lim - max_lim) > 1e-4
        lambda = (min_lim + max_lim)/2;
        fun_val = fi_fun(lambda, A, D, b, f);
        if  fun_val > 0
            min_lim = lambda;
        else
            max_lim = lambda;
        end
        num_iter = num_iter + 1;
    end
end

function fi = fi_fun(lambda_1, A, D, b, f)
    y = (A' * A + lambda_1 * D + 1e-6 * eye(size(D,1))) \ (A' * b - lambda_1 * f);
    fi = y' * D * y + 2 * f' * y;
end