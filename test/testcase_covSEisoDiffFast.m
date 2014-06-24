%% Setting
testcase_name = '[TEST_CASE: covSEisoDiffFast]';
disp(testcase_name);


%% Test: covSEisoDiffFast vs. covSEiso
test_name = '\t[TESE: covSEisoDiff vs. covSEiso]\n';
fprintf(1, test_name);

% K = K(x, x)
TEST_EQ(covSEisoDiffFast(hyp, x, [], 0, []), ...
        covSEiso(hyp, x), ...
        'K = K(x, x)');

% K_i = dK(x, x)/dtheta_i
for i = 1:length(hyp)
    TEST_EQ(covSEisoDiffFast(hyp, x, [], i, []), ...
            covSEiso(hyp, x, [], i), ...
            ['K_i = dK(x, x)/dtheta_', num2str(i)]);
end

% Ks = K(x, z)
TEST_EQ(covSEisoDiffFast(hyp, x, z, 0, []), ...
        covSEiso(hyp, x, z), ...
        'Ks = K(x, z)');

% Kss = K(z, z)
TEST_EQ(covSEisoDiffFast(hyp, z, 'diag', 0, []), ...
        covSEiso(hyp, z, 'diag'), ...
        'Kss = K(z, z)');


%% Test: covSEisoDiffFast vs. covSEisoDiffBwise
test_name = '\t[TESE: covSEisoDiffFast vs. covSEisoDiffBwise]\n';
fprintf(1, test_name);

% K = K(x, x)
TEST_EQ(covSEisoDiffFast(hyp, x, [], 0, xd), ...
        covSEisoDiff(hyp, x, [], 0, xd, true), ...
        'K = K(x, x)');

% K_i = dK(x, x)/dtheta_i    
for i = 1:length(hyp)                     
    TEST_EQ(covSEisoDiffFast(hyp, x, [], i, xd), ...
            covSEisoDiff(hyp, x, [], i, xd, true), ...
            ['K_i = dK(x, x)/dtheta_', num2str(i)]);
end

% Ks = K(x, z)
TEST_EQ(covSEisoDiffFast(hyp, x, z, 0, xd), ...
        covSEisoDiff(hyp, x, z, 0, xd, true), ...
        'Ks = K(x, z)');

% Kss = K(z, z)
TEST_EQ(covSEisoDiffFast(hyp, z, 'diag', 0, []), ...
        covSEisoDiff(hyp, z, 'diag', 0, [], true), ...
        'Kss = K(z, z)');

    
%% Test: covSEisoDiffFast - symmetric
test_name = '\t[TESE: covSEisoDiffFast - symmetric]\n';
fprintf(1, test_name);

% K = K(x, x)
K = covSEisoDiffFast(hyp, x, [], 0, xd);
TEST_EQ(K, ...
        K', ...
        'K = K(x, x)');

% K_i = dK(x, x)/dtheta_i    
for i = 1:length(hyp)
    K = covSEisoDiffFast(hyp, x, [], i, xd);
    TEST_EQ(K, ...
            K', ...
            ['K_i = dK(x, x)/dtheta_', num2str(i)]);
end

%% Test: covSEisoDiffFast - positive definite
test_name = '\t[TESE: covSEisoDiffFast - positive definite]\n';
fprintf(1, test_name);

% sigma_n
sigma_n = 1;

% K = K(x, x)
K = covSEisoDiffFast(hyp, x, [], 0, xd);
K = K + sigma_n^2*eye(size(K));
[R, p] = chol(K);
TEST_EQ(p, ...
        0', ...
        'K = K(x, x)');

% K_i = dK(x, x)/dtheta_i    
for i = 1:length(hyp)
    K = covSEisoDiffFast(hyp, x, [], i, xd);
    K = K + sigma_n^2*eye(size(K));
    [R, p] = chol(K);
    TEST_EQ(p, ...
            0', ...
            ['K_i = dK(x, x)/dtheta_', num2str(i)]);
end