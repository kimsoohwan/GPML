function testcase_covDerObs(f_covDerObs, f_covisoDerObsUnstable, f_cov, str_cov, hyp, x, z, xx, xxd, sigma_n)

% f_covDerObs:            function handle for the covariance function with derivative observations (stable version)
% f_covisoDerObsUnstable:    function handle for the covariance function with derivative observations (unstable version)
% f_cov:                function handle for the covariance function
% str_cov:              covariance function name
% hyp:                  hyperparameters
% x:                    first function inputs
% xd:                   first derivative inputs
% z:                    second function inputs
% sigma_n:              output noise variance

%% Setting
testcase_name = ['[TEST_CASE: ', str_cov, 'DerObs]'];
disp(testcase_name);


%% Test: cov*DerObs vs. cov*
test_name = ['\t[TESE: ', str_cov, 'DerObs vs. ', str_cov, ']\n'];
fprintf(1, test_name);

% K = K(x, x)
TEST_EQ(feval(f_covDerObs, hyp, xx, []), ...
        feval(f_cov{:}, hyp, x), ...
        'K = K(x, x)');

% K_i = dK(x, x)/dtheta_i
for i = 1:length(hyp)
    TEST_EQ(feval(f_covDerObs, hyp, xx, [], i), ...
            feval(f_cov{:}, hyp, x, [], i), ...
            ['K_', num2str(i), ' = dK(x, x)/dtheta_', num2str(i)]);
end

% Ks = K(x, z)
TEST_EQ(feval(f_covDerObs, hyp, xx, z), ...
        feval(f_cov{:}, hyp, x, z), ...
        'Ks = K(x, z)');

% Kss = K(z, z)
TEST_EQ(feval(f_covDerObs, hyp, z, 'diag'), ...
        feval(f_cov{:}, hyp, z, 'diag'), ...
        'Kss = K(z, z)');


%% Test: cov*DerObs vs. cov*DerObsBwiseUnstable
test_name = ['\t[TESE: ', str_cov, 'DerObs vs. ', str_cov, 'DerObsBwiseUnstable]\n'];
fprintf(1, test_name);

f_Bwise = true;

% K = K(x, x)
TEST_EQ(feval(f_covDerObs, hyp, xxd, [], 0), ...
        feval(f_covisoDerObsUnstable, hyp, xxd, [], 0, f_Bwise), ...
        'K = K(x, x)');

% K_i = dK(x, x)/dtheta_i    
for i = 1:length(hyp)                     
    TEST_EQ(feval(f_covDerObs, hyp, xxd, [], i), ...
            feval(f_covisoDerObsUnstable, hyp, xxd, [], i, f_Bwise), ...
            ['K_', num2str(i), ' = dK(x, x)/dtheta_', num2str(i)]);
end

% Ks = K(x, z)
TEST_EQ(feval(f_covDerObs, hyp, xxd, z, 0), ...
        feval(f_covisoDerObsUnstable, hyp, xxd, z, 0, f_Bwise), ...
        'Ks = K(x, z)');

% Kss = K(z, z)
TEST_EQ(feval(f_covDerObs, hyp, z, 'diag', 0), ...
        feval(f_covisoDerObsUnstable, hyp, z, 'diag', 0, f_Bwise), ...
        'Kss = K(z, z)');

    
%% Test: cov*DerObs - symmetric
test_name = ['\t[TESE: ', str_cov, 'DerObs - symmetric]\n'];
fprintf(1, test_name);

% K = K(x, x)
K = feval(f_covDerObs, hyp, xxd, [], 0);
TEST_EQ(K, ...
        K', ...
        'K = K(x, x)');

% K_i = dK(x, x)/dtheta_i    
for i = 1:length(hyp)
    K = feval(f_covDerObs, hyp, xxd, [], i);
    TEST_EQ(K, ...
            K', ...
            ['K_', num2str(i), ' = dK(x, x)/dtheta_', num2str(i)]);
end

%% Test: cov*DerObs - positive definite
test_name = ['\t[TESE: ', str_cov, 'DerObs - positive definite]\n'];
fprintf(1, test_name);

% K = K(x, x)
K = feval(f_covDerObs, hyp, xxd, [], 0);
K = K + sigma_n^2*eye(size(K));
[R, p] = chol(K);
TEST_EQ(p, ...
        0', ...
        'K = K(x, x)');

% % K_i = dK(x, x)/dtheta_i    
% for i = 1:length(hyp)
%     K = feval(f_covDerObs, hyp, xxd, [], i);
%     K = K + sigma_n^2*eye(size(K));
%     [R, p] = chol(K);
%     TEST_EQ(p, ...
%             0', ...
%             ['K_i = dK(x, x)/dtheta_', num2str(i)]);
% end

disp(' ');