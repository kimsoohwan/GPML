%% Setting
% hyperparameters
hyp = log([1, 1]);

% data
d = 3;
n  = 5;  x  = rand(n, d);
nd = 4;  xd = rand(nd, d);
ns = 3;  z  = rand(ns, d);

testcase_name = '[TEST_CASE: covMaternisoDiff]';
disp(testcase_name);


%% Test: covMaternisoDiffCwise vs. covMaterniso
test_name = '\t[TESE: covMaternisoDiffCwise vs. covMaterniso]\n';
fprintf(1, test_name);

% K = K(x, x)
TEST_EQ(covMaternisoDiff(hyp, x, [], 0, [], false), ...
        covMaterniso(3, hyp, x), ...
        'K = K(x, x)');

% K_i = dK(x, x)/dtheta_i
for i = 1:length(hyp)
    TEST_EQ(covMaternisoDiff(hyp, x, [], i, [], false), ...
            covMaterniso(3, hyp, x, [], i), ...
            ['K_i = dK(x, x)/dtheta_', num2str(i)]);
end

% Ks = K(x, z)
TEST_EQ(covMaternisoDiff(hyp, x, z, 0, [], false), ...
        covMaterniso(3, hyp, x, z), ...
        'Ks = K(x, z)');

% Kss = K(z, z)
TEST_EQ(covMaternisoDiff(hyp, z, 'diag', 0, [], false), ...
        covMaterniso(3, hyp, z, 'diag'), ...
        'Kss = K(z, z)');


%% Test: covMaternisoDiffBwise vs. covMaterniso
test_name = '\t[TESE: covMaternisoDiffBwise vs. covMaterniso]\n';
fprintf(1, test_name);

% K = K(x, x)
TEST_EQ(covMaternisoDiff(hyp, x, [], 0, [], true), ...
        covMaterniso(3, hyp, x), ...
        'K = K(x, x)');

% K_i = dK(x, x)/dtheta_i
for i = 1:length(hyp)
    TEST_EQ(covMaternisoDiff(hyp, x, [], i, [], true), ...
            covMaterniso(3, hyp, x, [], i), ...
            ['K_i = dK(x, x)/dtheta_', num2str(i)]);
end

% Ks = K(x, z)
TEST_EQ(covMaternisoDiff(hyp, x, z, 0, [], true), ...
        covMaterniso(3, hyp, x, z), ...
        'Ks = K(x, z)');

% Kss = K(z, z)
TEST_EQ(covMaternisoDiff(hyp, z, 'diag', 0, [], true), ...
        covMaterniso(3, hyp, z, 'diag'), ...
        'Kss = K(z, z)');


%% Test: covMaternisoDiffCwise vs. covMaternisoDiffBwise
test_name = '\t[TESE: covMaternisoDiffCwise vs. covMaternisoDiffBwise]\n';
fprintf(1, test_name);

% K = K(x, x)
TEST_EQ(covMaternisoDiff(hyp, x, [], 0, xd, false), ...
        covMaternisoDiff(hyp, x, [], 0, xd, true), ...
        'K = K(x, x)');

% K_i = dK(x, x)/dtheta_i    
for i = 1:length(hyp)                     
    TEST_EQ(covMaternisoDiff(hyp, x, [], i, xd, false), ...
            covMaternisoDiff(hyp, x, [], i, xd, true), ...
            ['K_i = dK(x, x)/dtheta_', num2str(i)]);
end

% Ks = K(x, z)
TEST_EQ(covMaternisoDiff(hyp, x, z, 0, xd, false), ...
        covMaternisoDiff(hyp, x, z, 0, xd, true), ...
        'Ks = K(x, z)');

% Kss = K(z, z)
TEST_EQ(covMaternisoDiff(hyp, z, 'diag', 0, [], false), ...
        covMaternisoDiff(hyp, z, 'diag', 0, [], true), ...
        'Kss = K(z, z)');

    
%% Test: covMaternisoDiffBwise - symmetric
test_name = '\t[TESE: covMaternisoDiffBwise - symmetric]\n';
fprintf(1, test_name);

% K = K(x, x)
K = covMaternisoDiff(hyp, x, [], 0, xd, true);
TEST_EQ(K, ...
        K', ...
        'K = K(x, x)');

% K_i = dK(x, x)/dtheta_i    
for i = 1:length(hyp)
    K = covMaternisoDiff(hyp, x, [], i, xd, true);
    TEST_EQ(K, ...
            K', ...
            ['K_i = dK(x, x)/dtheta_', num2str(i)]);
end

%% Test: covMaternisoDiffBwise - positive definite
test_name = '\t[TESE: covMaternisoDiffBwise - positive definite]\n';
fprintf(1, test_name);

% sigma_n
sigma_n = 1;

% K = K(x, x)
K = covMaternisoDiff(hyp, x, [], 0, xd, true);
K = K + sigma_n^2*eye(size(K));
[R, p] = chol(K);
TEST_EQ(p, ...
        0', ...
        'K = K(x, x)');

% K_i = dK(x, x)/dtheta_i    
for i = 1:length(hyp)
    K = covMaternisoDiff(hyp, x, [], i, xd, true);
    K = K + sigma_n^2*eye(size(K));
    [R, p] = chol(K);
    TEST_EQ(p, ...
            0', ...
            ['K_i = dK(x, x)/dtheta_', num2str(i)]);
end