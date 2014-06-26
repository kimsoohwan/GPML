function testcase_cov(f_cov, str_cov, hyp, x, sigma_n)

% f_cov:    function handle for the covariance function
% str_cov:  covariance function name
% hyp:                  hyperparameters
% x:                    first function inputs
% sigma_n:              output noise variance

%% Setting
testcase_name = ['[TEST_CASE: ', str_cov, ']'];
disp(testcase_name);

    
%% Test: cov* - symmetric
test_name = ['\t[TESE: ', str_cov, ' - symmetric]\n'];
fprintf(1, test_name);

% K = K(x, x)
K = feval(f_cov{:}, hyp, x);
TEST_EQ(K, ...
        K', ...
        'K = K(x, x)');

% K_i = dK(x, x)/dtheta_i    
for i = 1:length(hyp)
    K = feval(f_cov{:}, hyp, x, [], i);
    TEST_EQ(K, ...
            K', ...
            ['K_i = dK(x, x)/dtheta_', num2str(i)]);
end

%% Test: cov* - positive definite
test_name = ['\t[TESE: ', str_cov, ' - positive definite]\n'];
fprintf(1, test_name);

% K = K(x, x)
K = feval(f_cov{:}, hyp, x);
K = K + sigma_n^2*eye(size(K));
[R, p] = chol(K);
TEST_EQ(p, ...
        0', ...
        'K = K(x, x)');

% % K_i = dK(x, x)/dtheta_i    
% for i = 1:length(hyp)
%     K = feval(f_cov{:}, hyp, x, [], i);
%     K = K + sigma_n^2*eye(size(K));
%     [R, p] = chol(K);
%     TEST_EQ(p, ...
%             0', ...
%             ['K_i = dK(x, x)/dtheta_', num2str(i)]);
% end

disp(' ');