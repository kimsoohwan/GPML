%% Setting
testcase_name = '[TEST_CASE: covSparseisoDiffFast]';
disp(testcase_name);


%% Test: covSparseisoDiffFast vs. covSparseiso
test_name = '\t[TESE: covSparseisoDiff vs. covSparseiso]\n';
fprintf(1, test_name);

% K = K(x, x)
TEST_EQ(covSparseisoDiffFast(hyp, x, [], 0, []), ...
        covSparseiso(hyp, x), ...
        'K = K(x, x)');

% K_i = dK(x, x)/dtheta_i
for i = 1:length(hyp)
    TEST_EQ(covSparseisoDiffFast(hyp, x, [], i, []), ...
            covSparseiso(hyp, x, [], i), ...
            ['K_i = dK(x, x)/dtheta_', num2str(i)]);
end

% Ks = K(x, z)
TEST_EQ(covSparseisoDiffFast(hyp, x, z, 0, []), ...
        covSparseiso(hyp, x, z), ...
        'Ks = K(x, z)');

% Kss = K(z, z)
TEST_EQ(covSparseisoDiffFast(hyp, z, 'diag', 0, []), ...
        covSparseiso(hyp, z, 'diag'), ...
        'Kss = K(z, z)');


%% Test: covSparseisoDiffFast vs. covSparseisoDiffBwise
test_name = '\t[TESE: covSparseisoDiffFast vs. covSparseisoDiffBwise]\n';
fprintf(1, test_name);

% K = K(x, x)
TEST_EQ(covSparseisoDiffFast(hyp, x, [], 0, xd), ...
        covSparseisoDiff(hyp, x, [], 0, xd, true), ...
        'K = K(x, x)');

% K_i = dK(x, x)/dtheta_i    
for i = 1:length(hyp)                     
    TEST_EQ(covSparseisoDiffFast(hyp, x, [], i, xd), ...
            covSparseisoDiff(hyp, x, [], i, xd, true), ...
            ['K_i = dK(x, x)/dtheta_', num2str(i)]);
end

% Ks = K(x, z)
TEST_EQ(covSparseisoDiffFast(hyp, x, z, 0, xd), ...
        covSparseisoDiff(hyp, x, z, 0, xd, true), ...
        'Ks = K(x, z)');

% Kss = K(z, z)
TEST_EQ(covSparseisoDiffFast(hyp, z, 'diag', 0, []), ...
        covSparseisoDiff(hyp, z, 'diag', 0, [], true), ...
        'Kss = K(z, z)');

    
%% Test: covSparseisoDiffFast - symmetric
test_name = '\t[TESE: covSparseisoDiffFast - symmetric]\n';
fprintf(1, test_name);

% K = K(x, x)
K = covSparseisoDiffFast(hyp, x, [], 0, xd);
TEST_EQ(K, ...
        K', ...
        'K = K(x, x)');

% K_i = dK(x, x)/dtheta_i    
for i = 1:length(hyp)
    K = covSparseisoDiffFast(hyp, x, [], i, xd);
    TEST_EQ(K, ...
            K', ...
            ['K_i = dK(x, x)/dtheta_', num2str(i)]);
end

%% Test: covSparseisoDiffFast - positive definite
test_name = '\t[TESE: covSparseisoDiffFast - positive definite]\n';
fprintf(1, test_name);

% sigma_n
sigma_n = 1;

% K = K(x, x)
K = covSparseisoDiffFast(hyp, x, [], 0, xd);
K = K + sigma_n^2*eye(size(K));
[R, p] = chol(K);
TEST_EQ(p, ...
        0', ...
        'K = K(x, x)');

% K_i = dK(x, x)/dtheta_i    
for i = 1:length(hyp)
    K = covSparseisoDiffFast(hyp, x, [], i, xd);
    K = K + sigma_n^2*eye(size(K));
    [R, p] = chol(K);
    TEST_EQ(p, ...
            0', ...
            ['K_i = dK(x, x)/dtheta_', num2str(i)]);
end