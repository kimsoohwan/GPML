function testcase_covisoDerObsUnstable(f_covisoDerObsUnstable, f_cov, str_cov, hyp, x, z, xx, xxd, sigma_n)

% f_covisoDerObsUnstable:    function handle for the covariance function with derivative observations (unstable version)
% f_cov:                function handle for the covariance function
% str_cov:              covariance function name
% hyp:                  hyperparameters
% x:                    first function inputs
% xd:                   first derivative inputs
% z:                    second function inputs
% sigma_n:              output noise variance

%% Setting
testcase_name = ['[TEST_CASE: ', str_cov, 'DerObsUnstable]'];
disp(testcase_name);

f_Cwise = false;
f_Bwise = true;


%% Test: cov*DerObsCwiseUnstable vs. cov*
test_name = ['\t[TESE: ', str_cov, 'DerObsCwiseUnstable vs. ', str_cov, ']\n'];
fprintf(1, test_name);

% K = K(x, x)
TEST_EQ(feval(f_covisoDerObsUnstable, hyp, xx, [], 0, f_Cwise), ...
        feval(f_cov{:}, hyp, x), ...
        'K = K(x, x)');

% K_i = dK(x, x)/dtheta_i
for i = 1:length(hyp)
    TEST_EQ(feval(f_covisoDerObsUnstable, hyp, xx, [], i, f_Cwise), ...
            feval(f_cov{:}, hyp, x, [], i), ...
            ['K_', num2str(i), ' = dK(x, x)/dtheta_', num2str(i)]);
end

% Ks = K(x, z)
TEST_EQ(feval(f_covisoDerObsUnstable, hyp, xx, z, 0, f_Cwise), ...
        feval(f_cov{:}, hyp, x, z), ...
        'Ks = K(x, z)');

% Kss = K(z, z)
TEST_EQ(feval(f_covisoDerObsUnstable, hyp, z, 'diag', 0, f_Cwise), ...
        feval(f_cov{:}, hyp, z, 'diag'), ...
        'Kss = K(z, z)');


%% Test: cov*DerObsBwiseUnstable vs. cov*
test_name = ['\t[TESE: ', str_cov, 'DerObsBwiseUnstable vs. ', str_cov, ']\n'];
fprintf(1, test_name);

% K = K(x, x)
TEST_EQ(feval(f_covisoDerObsUnstable, hyp, xx, [], 0, f_Bwise), ...
        feval(f_cov{:}, hyp, x), ...
        'K = K(x, x)');

% K_i = dK(x, x)/dtheta_i
for i = 1:length(hyp)
    TEST_EQ(feval(f_covisoDerObsUnstable, hyp, xx, [], i, f_Bwise), ...
            feval(f_cov{:}, hyp, x, [], i), ...
            ['K_', num2str(i), ' = dK(x, x)/dtheta_', num2str(i)]);
end

% Ks = K(x, z)
TEST_EQ(feval(f_covisoDerObsUnstable, hyp, xx, z, 0, f_Bwise), ...
        feval(f_cov{:}, hyp, x, z), ...
        'Ks = K(x, z)');

% Kss = K(z, z)
TEST_EQ(feval(f_covisoDerObsUnstable, hyp, z, 'diag', 0, f_Bwise), ...
        feval(f_cov{:}, hyp, z, 'diag'), ...
        'Kss = K(z, z)');


%% Test: cov*DerObsCwiseUnstable vs. cov*DerObsBwiseUnstable
test_name = ['\t[TESE: ', str_cov, 'DerObsCwiseUnstable vs. ', str_cov, 'DerObsBwiseUnstable]\n'];
fprintf(1, test_name);

% K = K(x, x)
TEST_EQ(feval(f_covisoDerObsUnstable, hyp, xxd, [], 0, f_Cwise), ...
        feval(f_covisoDerObsUnstable, hyp, xxd, [], 0, f_Bwise), ...
        'K = K(x, x)');

% K_i = dK(x, x)/dtheta_i    
for i = 1:length(hyp)                     
    TEST_EQ(feval(f_covisoDerObsUnstable, hyp, xxd, [], i, f_Cwise), ...
            feval(f_covisoDerObsUnstable, hyp, xxd, [], i, f_Bwise), ...
            ['K_', num2str(i), ' = dK(x, x)/dtheta_', num2str(i)]);
end

% Ks = K(x, z)
TEST_EQ(feval(f_covisoDerObsUnstable, hyp, xxd, z, 0, f_Cwise), ...
        feval(f_covisoDerObsUnstable, hyp, xxd, z, 0, f_Bwise), ...
        'Ks = K(x, z)');

% Kss = K(z, z)
TEST_EQ(feval(f_covisoDerObsUnstable, hyp, z, 'diag', 0, f_Cwise), ...
        feval(f_covisoDerObsUnstable, hyp, z, 'diag', 0, f_Bwise), ...
        'Kss = K(z, z)');

    
%% Test: cov*DerObsBwiseUnstable - symmetric
test_name = ['\t[TESE: ', str_cov, 'DerObsBwiseUnstable - symmetric]\n'];
fprintf(1, test_name);

% K = K(x, x)
K = feval(f_covisoDerObsUnstable, hyp, xxd, [], 0, f_Bwise);
TEST_EQ(K, ...
        K', ...
        'K = K(x, x)');

% K_i = dK(x, x)/dtheta_i    
for i = 1:length(hyp)
    K = feval(f_covisoDerObsUnstable, hyp, xxd, [], i, f_Bwise);
    TEST_EQ(K, ...
            K', ...
            ['K_', num2str(i), ' = dK(x, x)/dtheta_', num2str(i)]);
end

%% Test: cov*DerObsBwiseUnstable - positive definite
test_name = ['\t[TESE: ', str_cov, 'DerObsBwiseUnstable - positive definite]\n'];
fprintf(1, test_name);

% K = K(x, x)
K = feval(f_covisoDerObsUnstable, hyp, xxd, [], 0, f_Bwise);
K = K + sigma_n^2*eye(size(K));
[R, p] = chol(K);
TEST_EQ(p, ...
        0', ...
        'K = K(x, x)');

% % K_i = dK(x, x)/dtheta_i    
% for i = 1:length(hyp)
%     K = feval(f_covisoDerObsUnstable, hyp, xxd, [], i, f_Bwise);
%     K = K + sigma_n^2*eye(size(K));
%     [R, p] = chol(K);
%     TEST_EQ(p, ...
%             0', ...
%             ['K_i = dK(x, x)/dtheta_', num2str(i)]);
% end

disp(' ');