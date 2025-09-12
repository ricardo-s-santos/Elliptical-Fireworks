% Ovo je kod za fi funkciju iz Beck-ovog rada koju koristim za bisekciju

function fi = fi_fun(lambda_1, A, D, b, f)

% if rcond(inv(A' * A + lambda_1 * D)) > 10^(-4)

%     y = (A' * A + lambda_1 * D) \ (A' * b - lambda_1 * f);
    y = (A' * A + lambda_1 * D + 1e-6 * eye(size(D,1))) \ (A' * b - lambda_1 * f);

    fi = y' * D * y + 2 * f' * y;
    
% end