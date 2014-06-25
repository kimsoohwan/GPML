%% Setting
testcase_name = '[TEST_CASE: covSparseisoDiff]';
disp(testcase_name);


%% Test: covSparseisoDiffCwise vs. covSparseiso
test_name = '\t[TESE: covSparseisoDiffCwise vs. covSparseiso]\n';
fprintf(1, test_name);

% K = K(x, x)
TEST_EQ(covSparseisoDiff(hyp, x, [], 0, [], false), ...
        covSparseiso(hyp, x), ...
        'K = K(x, x)');

% K_i = dK(x, x)/dtheta_i
for i = 1:length(hyp)
    TEST_EQ(covSparseisoDiff(hyp, x, [], i, [], false), ...
            covSparseiso(hyp, x, [], i), ...
            ['K_i = dK(x, x)/dtheta_', num2str(i)]);
end

% Ks = K(x, z)
TEST_EQ(covSparseisoDiff(hyp, x, z, 0, [], false), ...
        covSparseiso(hyp, x, z), ...
        'Ks = K(x, z)');

% Kss = K(z, z)
TEST_EQ(covSparseisoDiff(hyp, z, 'diag', 0, [], false), ...
        covSparseiso(hyp, z, 'diag'), ...
        'Kss = K(z, z)');


%% Test: covSparseisoDiffBwise vs. covSparseiso
test_name = '\t[TESE: covSparseisoDiffBwise vs. covSparseiso]\n';
fprintf(1, test_name);

% K = K(x, x)
TEST_EQ(covSparseisoDiff(hyp, x, [], 0, [], true), ...
        covSparseiso(hyp, x), ...
        'K = K(x, x)');

% K_i = dK(x, x)/dtheta_i
for i = 1:length(hyp)
    TEST_EQ(covSparseisoDiff(hyp, x, [], i, [], true), ...
            covSparseiso(hyp, x, [], i), ...
            ['K_i = dK(x, x)/dtheta_', num2str(i)]);
end

% Ks = K(x, z)
TEST_EQ(covSparseisoDiff(hyp, x, z, 0, [], true), ...
        covSparseiso(hyp, x, z), ...
        'Ks = K(x, z)');

% Kss = K(z, z)
TEST_EQ(covSparseisoDiff(hyp, z, 'diag', 0, [], true), ...
        covSparseiso(hyp, z, 'diag'), ...
        'Kss = K(z, z)');


%% Test: covSparseisoDiffCwise vs. covSparseisoDiffBwise
test_name = '\t[TESE: covSparseisoDiffCwise vs. covSparseisoDiffBwise]\n';
fprintf(1, test_name);

% K = K(x, x)
TEST_EQ(covSparseisoDiff(hyp, x, [], 0, xd, false), ...
        covSparseisoDiff(hyp, x, [], 0, xd, true), ...
        'K = K(x, x)');

% K_i = dK(x, x)/dtheta_i    
for i = 1:length(hyp)                     
    TEST_EQ(covSparseisoDiff(hyp, x, [], i, xd, false), ...
            covSparseisoDiff(hyp, x, [], i, xd, true), ...
            ['K_i = dK(x, x)/dtheta_', num2str(i)]);
end

% Ks = K(x, z)
TEST_EQ(covSparseisoDiff(hyp, x, z, 0, xd, false), ...
        covSparseisoDiff(hyp, x, z, 0, xd, true), ...
        'Ks = K(x, z)');

% Kss = K(z, z)
TEST_EQ(covSparseisoDiff(hyp, z, 'diag', 0, [], false), ...
        covSparseisoDiff(hyp, z, 'diag', 0, [], true), ...
        'Kss = K(z, z)');

    
%% Test: covSparseisoDiffBwise - symmetric
test_name = '\t[TESE: covSparseisoDiffBwise - symmetric]\n';
fprintf(1, test_name);

% K = K(x, x)
K = covSparseisoDiff(hyp, x, [], 0, xd, true);
TEST_EQ(K, ...
        K', ...
        'K = K(x, x)');

% K_i = dK(x, x)/dtheta_i    
for i = 1:length(hyp)
    K = covSparseisoDiff(hyp, x, [], i, xd, true);
    TEST_EQ(K, ...
            K', ...
            ['K_i = dK(x, x)/dtheta_', num2str(i)]);
end

%% Test: covSparseisoDiffBwise - positive definite
test_name = '\t[TESE: covSparseisoDiffBwise - positive definite]\n';
fprintf(1, test_name);

% sigma_n
sigma_n = 1;

% K = K(x, x)
K = covSparseisoDiff(hyp, x, [], 0, xd, true);
K = K + sigma_n^2*eye(size(K));
[R, p] = chol(K);
TEST_EQ(p, ...
        0', ...
        'K = K(x, x)');

% K_i = dK(x, x)/dtheta_i    
for i = 1:length(hyp)
    K = covSparseisoDiff(hyp, x, [], i, xd, true);
    K = K + sigma_n^2*eye(size(K));
    [R, p] = chol(K);
    TEST_EQ(p, ...
            0', ...
            ['K_i = dK(x, x)/dtheta_', num2str(i)]);
end