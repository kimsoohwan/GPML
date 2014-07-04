function print_covDerObs(f_covDerObs, str_cov, hyp, xxd, z)

% f_covDerObs:            function handle for the covariance function with derivative observations (stable version)
% str_cov:              covariance function name
% hyp:                  hyperparameters
% x:                    first function inputs
% xd:                   first derivative inputs
% z:                    second function inputs

% Setting
printcase_name = ['[PRINT: ', str_cov, 'DerObs]'];
disp(printcase_name);

% K = K(x, x)
K = feval(f_covDerObs, hyp, xxd);
print_matrix_for_reference(K, 'K = K(x, x)');

% K_i = dK(x, x)/dtheta_i    
for i = 1:length(hyp)   
    K = feval(f_covDerObs, hyp, xxd, [], i);
    print_matrix_for_reference(K, ['K_', num2str(i), ' = dK(x, x)/dtheta_', num2str(i)]);
end

% Ks = K(x, z)
K = feval(f_covDerObs, hyp, xxd, z);
print_matrix_for_reference(K, 'Ks = K(x, z)');

% Kss = K(z, z)
K = feval(f_covDerObs, hyp, z, 'diag');
print_matrix_for_reference(K, 'Kss = K(z, z)');

disp(' ');