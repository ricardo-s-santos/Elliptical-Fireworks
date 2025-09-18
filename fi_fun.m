function fi = fi_fun(lambda_1, A, D, b, f)
    y = (A' * A + lambda_1 * D + 1e-6 * eye(size(D,1))) \ (A' * b - lambda_1 * f);
    fi = y' * D * y + 2 * f' * y;
end