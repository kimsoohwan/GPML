function testcase_covDiff(f_covDiff, f_covDiffUnstable, f_cov, str_cov, hyp, x, xd, z, sigma_n)

% f_covDiff:            function handle for the covariance function with derivative observations (stable version)
% f_covDiffUnstable:    function handle for the covariance function with derivative observations (unstable version)
% f_cov:                function handle for the covariance function
% str_cov:              covariance function name
% hyp:                  hyperparameters
% x:                    first function inputs
% xd:                   first derivative inputs
% z:                    second function inputs
% sigma_n:              output noise variance

%% Setting
testcase_name = ['[TEST_CASE: ', str_cov, 'Diff]'];
disp(testcase_name);


%% Test: cov*Diff vs. cov*
test_name = ['\t[TESE: ', str_cov, 'Diff vs. ', str_cov, ']\n'];
fprintf(1, test_name);

% K = K(x, x)
TEST_EQ(feval(f_covDiff, hyp, x, [], 0, []), ...
        feval(f_cov{:}, hyp, x), ...
        'K = K(x, x)');

% K_i = dK(x, x)/dtheta_i
for i = 1:length(hyp)
    TEST_EQ(feval(f_covDiff, hyp, x, [], i, []), ...
            feval(f_cov{:}, hyp, x, [], i), ...
            ['K_i = dK(x, x)/dtheta_', num2str(i)]);
end

% Ks = K(x, z)
TEST_EQ(feval(f_covDiff, hyp, x, z, 0, []), ...
        feval(f_cov{:}, hyp, x, z), ...
        'Ks = K(x, z)');

% Kss = K(z, z)
TEST_EQ(feval(f_covDiff, hyp, z, 'diag', 0, []), ...
        feval(f_cov{:}, hyp, z, 'diag'), ...
        'Kss = K(z, z)');


%% Test: cov*Diff vs. cov*DiffBwiseUnstable
test_name = ['\t[TESE: ', str_cov, 'Diff vs. ', str_cov, 'DiffBwiseUnstable]\n'];
fprintf(1, test_name);

f_Bwise = true;

% K = K(x, x)
TEST_EQ(feval(f_covDiff, hyp, x, [], 0, xd), ...
        feval(f_covDiffUnstable, hyp, x, [], 0, xd, f_Bwise), ...
        'K = K(x, x)');

% K_i = dK(x, x)/dtheta_i    
for i = 1:length(hyp)                     
    TEST_EQ(feval(f_covDiff, hyp, x, [], i, xd), ...
            feval(f_covDiffUnstable, hyp, x, [], i, xd, f_Bwise), ...
            ['K_i = dK(x, x)/dtheta_', num2str(i)]);
end

% Ks = K(x, z)
TEST_EQ(feval(f_covDiff, hyp, x, z, 0, xd), ...
        feval(f_covDiffUnstable, hyp, x, z, 0, xd, f_Bwise), ...
        'Ks = K(x, z)');

% Kss = K(z, z)
TEST_EQ(feval(f_covDiff, hyp, z, 'diag', 0, []), ...
        feval(f_covDiffUnstable, hyp, z, 'diag', 0, [], f_Bwise), ...
        'Kss = K(z, z)');

    
%% Test: cov*Diff - symmetric
test_name = ['\t[TESE: ', str_cov, 'Diff - symmetric]\n'];
fprintf(1, test_name);

% K = K(x, x)
K = feval(f_covDiff, hyp, x, [], 0, xd);
TEST_EQ(K, ...
        K', ...
        'K = K(x, x)');

% K_i = dK(x, x)/dtheta_i    
for i = 1:length(hyp)
    K = feval(f_covDiff, hyp, x, [], i, xd);
    TEST_EQ(K, ...
            K', ...
            ['K_i = dK(x, x)/dtheta_', num2str(i)]);
end

%% Test: cov*Diff - positive definite
test_name = ['\t[TESE: ', str_cov, 'Diff - positive definite]\n'];
fprintf(1, test_name);

% K = K(x, x)
K = feval(f_covDiff, hyp, x, [], 0, xd);
K = K + sigma_n^2*eye(size(K));
[R, p] = chol(K);
TEST_EQ(p, ...
        0', ...
        'K = K(x, x)');

% % K_i = dK(x, x)/dtheta_i    
% for i = 1:length(hyp)
%     K = feval(f_covDiff, hyp, x, [], i, xd);
%     K = K + sigma_n^2*eye(size(K));
%     [R, p] = chol(K);
%     TEST_EQ(p, ...
%             0', ...
%             ['K_i = dK(x, x)/dtheta_', num2str(i)]);
% end

disp(' ');