%% Setting
testcase_name = '[TEST_CASE: covMaterniso3DiffFast]';
disp(testcase_name);


%% Test: covMaterniso3DiffFast vs. covMaterniso
test_name = '\t[TESE: covMaterniso3DiffFast vs. covMaterniso]\n';
fprintf(1, test_name);

% K = K(x, x)
TEST_EQ(covMaterniso3DiffFast(hyp, x, [], 0, []), ...
        covMaterniso(3, hyp, x), ...
        'K = K(x, x)');

% K_i = dK(x, x)/dtheta_i
for i = 1:length(hyp)
    TEST_EQ(covMaterniso3DiffFast(hyp, x, [], i, []), ...
            covMaterniso(3, hyp, x, [], i), ...
            ['K_i = dK(x, x)/dtheta_', num2str(i)]);
end

% Ks = K(x, z)
TEST_EQ(covMaterniso3DiffFast(hyp, x, z, 0, []), ...
        covMaterniso(3, hyp, x, z), ...
        'Ks = K(x, z)');

% Kss = K(z, z)
TEST_EQ(covMaterniso3DiffFast(hyp, z, 'diag', 0, []), ...
        covMaterniso(3, hyp, z, 'diag'), ...
        'Kss = K(z, z)');


%% Test: covMaterniso3DiffFast vs. covMaterniso3DiffBwise
test_name = '\t[TESE: covMaterniso3DiffFast vs. covMaterniso3DiffBwise]\n';
fprintf(1, test_name);

% K = K(x, x)
TEST_EQ(covMaterniso3DiffFast(hyp, x, [], 0, xd), ...
        covMaterniso3Diff(hyp, x, [], 0, xd, true), ...
        'K = K(x, x)');

% K_i = dK(x, x)/dtheta_i    
for i = 1:length(hyp)                     
    TEST_EQ(covMaterniso3DiffFast(hyp, x, [], i, xd), ...
            covMaterniso3Diff(hyp, x, [], i, xd, true), ...
            ['K_i = dK(x, x)/dtheta_', num2str(i)]);
end

% Ks = K(x, z)
TEST_EQ(covMaterniso3DiffFast(hyp, x, z, 0, xd), ...
        covMaterniso3Diff(hyp, x, z, 0, xd, true), ...
        'Ks = K(x, z)');

% Kss = K(z, z)
TEST_EQ(covMaterniso3DiffFast(hyp, z, 'diag', 0, []), ...
        covMaterniso3Diff(hyp, z, 'diag', 0, [], true), ...
        'Kss = K(z, z)');

    
%% Test: covMaterniso3DiffBwise - symmetric
test_name = '\t[TESE: covMaterniso3DiffFast - symmetric]\n';
fprintf(1, test_name);

% K = K(x, x)
K = covMaterniso3DiffFast(hyp, x, [], 0, xd);
TEST_EQ(K, ...
        K', ...
        'K = K(x, x)');

% K_i = dK(x, x)/dtheta_i    
for i = 1:length(hyp)
    K = covMaterniso3DiffFast(hyp, x, [], i, xd);
    TEST_EQ(K, ...
            K', ...
            ['K_i = dK(x, x)/dtheta_', num2str(i)]);
end

%% Test: covMaterniso3DiffFast - positive definite
test_name = '\t[TESE: covMaterniso3DiffFast - positive definite]\n';
fprintf(1, test_name);

% sigma_n
sigma_n = 1;

% K = K(x, x)
K = covMaterniso3DiffFast(hyp, x, [], 0, xd);
K = K + sigma_n^2*eye(size(K));
[R, p] = chol(K);
TEST_EQ(p, ...
        0', ...
        'K = K(x, x)');

% K_i = dK(x, x)/dtheta_i    
for i = 1:length(hyp)
    K = covMaterniso3DiffFast(hyp, x, [], i, xd);
    K = K + sigma_n^2*eye(size(K));
    [R, p] = chol(K);
    TEST_EQ(p, ...
            0', ...
            ['K_i = dK(x, x)/dtheta_', num2str(i)]);
end