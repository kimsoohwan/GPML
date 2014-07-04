function print_cov(f_cov, str_cov, hyp, x, z)

% f_cov:                function handle for the covariance function
% str_cov:              covariance function name
% hyp:                  hyperparameters
% x:                    first function inputs
% z:                    second function inputs

% Setting
printcase_name = ['[PRINT: ', str_cov, ']'];
disp(printcase_name);

% K = K(x, x)
K = feval(f_cov{:}, hyp, x);
print_matrix_for_reference(K, 'K = K(x, x)');

% K_i = dK(x, x)/dtheta_i
for i = 1:length(hyp)
    K = feval(f_cov{:}, hyp, x, [], i);
    print_matrix_for_reference(K, ['K_', num2str(i), ' = dK(x, x)/dtheta_', num2str(i)]);
end

% Ks = K(x, z)
K = feval(f_cov{:}, hyp, x, z);
print_matrix_for_reference(K, 'Ks = K(x, z)');

% Kss = K(z, z)
K = feval(f_cov{:}, hyp, z, 'diag');
print_matrix_for_reference(K, 'Kss = K(z, z)');

disp(' ');