%% Setting
testcase_name = '[TEST_CASE: covMaterniso3Diff]';
disp(testcase_name);


%% Test: covMaterniso3DiffCwise vs. covMaterniso
test_name = '\t[TESE: covMaterniso3DiffCwise vs. covMaterniso]\n';
fprintf(1, test_name);

% K = K(x, x)
TEST_EQ(covMaterniso3Diff(hyp, x, [], 0, [], false), ...
        covMaterniso(3, hyp, x), ...
        'K = K(x, x)');

% K_i = dK(x, x)/dtheta_i
for i = 1:length(hyp)
    TEST_EQ(covMaterniso3Diff(hyp, x, [], i, [], false), ...
            covMaterniso(3, hyp, x, [], i), ...
            ['K_i = dK(x, x)/dtheta_', num2str(i)]);
end

% Ks = K(x, z)
TEST_EQ(covMaterniso3Diff(hyp, x, z, 0, [], false), ...
        covMaterniso(3, hyp, x, z), ...
        'Ks = K(x, z)');

% Kss = K(z, z)
TEST_EQ(covMaterniso3Diff(hyp, z, 'diag', 0, [], false), ...
        covMaterniso(3, hyp, z, 'diag'), ...
        'Kss = K(z, z)');


%% Test: covMaterniso3DiffBwise vs. covMaterniso
test_name = '\t[TESE: covMaterniso3DiffBwise vs. covMaterniso]\n';
fprintf(1, test_name);

% K = K(x, x)
TEST_EQ(covMaterniso3Diff(hyp, x, [], 0, [], true), ...
        covMaterniso(3, hyp, x), ...
        'K = K(x, x)');

% K_i = dK(x, x)/dtheta_i
for i = 1:length(hyp)
    TEST_EQ(covMaterniso3Diff(hyp, x, [], i, [], true), ...
            covMaterniso(3, hyp, x, [], i), ...
            ['K_i = dK(x, x)/dtheta_', num2str(i)]);
end

% Ks = K(x, z)
TEST_EQ(covMaterniso3Diff(hyp, x, z, 0, [], true), ...
        covMaterniso(3, hyp, x, z), ...
        'Ks = K(x, z)');

% Kss = K(z, z)
TEST_EQ(covMaterniso3Diff(hyp, z, 'diag', 0, [], true), ...
        covMaterniso(3, hyp, z, 'diag'), ...
        'Kss = K(z, z)');


%% Test: covMaterniso3DiffCwise vs. covMaterniso3DiffBwise
test_name = '\t[TESE: covMaterniso3DiffCwise vs. covMaterniso3DiffBwise]\n';
fprintf(1, test_name);

% K = K(x, x)
TEST_EQ(covMaterniso3Diff(hyp, x, [], 0, xd, false), ...
        covMaterniso3Diff(hyp, x, [], 0, xd, true), ...
        'K = K(x, x)');

% K_i = dK(x, x)/dtheta_i    
for i = 1:length(hyp)                     
    TEST_EQ(covMaterniso3Diff(hyp, x, [], i, xd, false), ...
            covMaterniso3Diff(hyp, x, [], i, xd, true), ...
            ['K_i = dK(x, x)/dtheta_', num2str(i)]);
end

% Ks = K(x, z)
TEST_EQ(covMaterniso3Diff(hyp, x, z, 0, xd, false), ...
        covMaterniso3Diff(hyp, x, z, 0, xd, true), ...
        'Ks = K(x, z)');

% Kss = K(z, z)
TEST_EQ(covMaterniso3Diff(hyp, z, 'diag', 0, [], false), ...
        covMaterniso3Diff(hyp, z, 'diag', 0, [], true), ...
        'Kss = K(z, z)');

    
%% Test: covMaterniso3DiffBwise - symmetric
test_name = '\t[TESE: covMaterniso3DiffBwise - symmetric]\n';
fprintf(1, test_name);

% K = K(x, x)
K = covMaterniso3Diff(hyp, x, [], 0, xd, true);
TEST_EQ(K, ...
        K', ...
        'K = K(x, x)');

% K_i = dK(x, x)/dtheta_i    
for i = 1:length(hyp)
    K = covMaterniso3Diff(hyp, x, [], i, xd, true);
    TEST_EQ(K, ...
            K', ...
            ['K_i = dK(x, x)/dtheta_', num2str(i)]);
end

%% Test: covMaterniso3DiffBwise - positive definite
test_name = '\t[TESE: covMaterniso3DiffBwise - positive definite]\n';
fprintf(1, test_name);

% sigma_n
sigma_n = 1;

% K = K(x, x)
K = covMaterniso3Diff(hyp, x, [], 0, xd, true);
K = K + sigma_n^2*eye(size(K));
[R, p] = chol(K);
TEST_EQ(p, ...
        0', ...
        'K = K(x, x)');

% K_i = dK(x, x)/dtheta_i    
for i = 1:length(hyp)
    K = covMaterniso3Diff(hyp, x, [], i, xd, true);
    K = K + sigma_n^2*eye(size(K));
    [R, p] = chol(K);
    TEST_EQ(p, ...
            0', ...
            ['K_i = dK(x, x)/dtheta_', num2str(i)]);
end