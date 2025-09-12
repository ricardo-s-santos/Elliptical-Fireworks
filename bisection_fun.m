% Ovo je kod za bisekciju koju koristim u Wang-ovom radu

function lambda = bisection_fun(min_lim, max_lim, tol, N_iter, A, D, b, f)

% if fi_fun(min_lim, A, D, b, f) <= tol % proveravam da li je vrednost funkcije u pocetnim tackama manja od tolerancije
%     
%     lamda = min_lim;
%     
%     return;
%     
% elseif fi_fun(max_lim, A, D, b, f) <= tol
%     
%     lamda = max_lim;
%     
%     return;
%     
% end


% Sada radim bisekciju

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