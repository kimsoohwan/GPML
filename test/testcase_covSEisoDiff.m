%% Setting
% hyperparameters
hyp = log([1, 1]);

% data
d = 3;
n  = 5;  x  = rand(n, d);
nd = 4;  xd = rand(nd, d);
ns = 3;  z  = rand(ns, d);

testcase_name = '[TEST_CASE: covSEisoDiff]';
disp(testcase_name);


%% Test: covSEisoDiffCwise vs. covSEiso
test_name = '\t[TESE: covSEisoDiffCwise vs. covSEiso]\n';
fprintf(1, test_name);

% K = K(x, x)
TEST_EQ(covSEisoDiff(hyp, x, [], 0, [], false), ...
        covSEiso(hyp, x), ...
        'K = K(x, x)');

% K_i = dK(x, x)/dtheta_i
for i = 1:length(hyp)
    TEST_EQ(covSEisoDiff(hyp, x, [], i, [], false), ...
            covSEiso(hyp, x, [], i), ...
            ['K_i = dK(x, x)/dtheta_', num2str(i)]);
end

% Ks = K(x, z)
TEST_EQ(covSEisoDiff(hyp, x, z, 0, [], false), ...
        covSEiso(hyp, x, z), ...
        'Ks = K(x, z)');

% Kss = K(z, z)
TEST_EQ(covSEisoDiff(hyp, z, 'diag', 0, [], false), ...
        covSEiso(hyp, z, 'diag'), ...
        'Kss = K(z, z)');


%% Test: covSEisoDiffBwise vs. covSEiso
test_name = '\t[TESE: covSEisoDiffBwise vs. covSEiso]\n';
fprintf(1, test_name);

% K = K(x, x)
TEST_EQ(covSEisoDiff(hyp, x, [], 0, [], true), ...
        covSEiso(hyp, x), ...
        'K = K(x, x)');

% K_i = dK(x, x)/dtheta_i
for i = 1:length(hyp)
    TEST_EQ(covSEisoDiff(hyp, x, [], i, [], true), ...
            covSEiso(hyp, x, [], i), ...
            ['K_i = dK(x, x)/dtheta_', num2str(i)]);
end

% Ks = K(x, z)
TEST_EQ(covSEisoDiff(hyp, x, z, 0, [], true), ...
        covSEiso(hyp, x, z), ...
        'Ks = K(x, z)');

% Kss = K(z, z)
TEST_EQ(covSEisoDiff(hyp, z, 'diag', 0, [], true), ...
        covSEiso(hyp, z, 'diag'), ...
        'Kss = K(z, z)');


%% Test: covSEisoDiffCwise vs. covSEisoDiffBwise
test_name = '\t[TESE: covSEisoDiffCwise vs. covSEisoDiffBwise]\n';
fprintf(1, test_name);

% K = K(x, x)
TEST_EQ(covSEisoDiff(hyp, x, [], 0, xd, false), ...
        covSEisoDiff(hyp, x, [], 0, xd, true), ...
        'K = K(x, x)');

% K_i = dK(x, x)/dtheta_i    
for i = 1:length(hyp)                     
    TEST_EQ(covSEisoDiff(hyp, x, [], i, xd, false), ...
            covSEisoDiff(hyp, x, [], i, xd, true), ...
            ['K_i = dK(x, x)/dtheta_', num2str(i)]);
end

% Ks = K(x, z)
TEST_EQ(covSEisoDiff(hyp, x, z, 0, xd, false), ...
        covSEisoDiff(hyp, x, z, 0, xd, true), ...
        'Ks = K(x, z)');

% Kss = K(z, z)
TEST_EQ(covSEisoDiff(hyp, z, 'diag', 0, [], false), ...
        covSEisoDiff(hyp, z, 'diag', 0, [], true), ...
        'Kss = K(z, z)');

    
%% Test: covSEisoDiffBwise - symmetric
test_name = '\t[TESE: covSEisoDiffBwise - symmetric]\n';
fprintf(1, test_name);

% K = K(x, x)
K = covSEisoDiff(hyp, x, [], 0, xd, true);
TEST_EQ(K, ...
        K', ...
        'K = K(x, x)');

% K_i = dK(x, x)/dtheta_i    
for i = 1:length(hyp)
    K = covSEisoDiff(hyp, x, [], i, xd, true);
    TEST_EQ(K, ...
            K', ...
            ['K_i = dK(x, x)/dtheta_', num2str(i)]);
end

%% Test: covSEisoDiffBwise - positive definite
test_name = '\t[TESE: covSEisoDiffBwise - positive definite]\n';
fprintf(1, test_name);

% sigma_n
sigma_n = 1;

% K = K(x, x)
K = covSEisoDiff(hyp, x, [], 0, xd, true);
K = K + sigma_n^2*eye(size(K));
[R, p] = chol(K);
TEST_EQ(p, ...
        0', ...
        'K = K(x, x)');

% K_i = dK(x, x)/dtheta_i    
for i = 1:length(hyp)
    K = covSEisoDiff(hyp, x, [], i, xd, true);
    K = K + sigma_n^2*eye(size(K));
    [R, p] = chol(K);
    TEST_EQ(p, ...
            0', ...
            ['K_i = dK(x, x)/dtheta_', num2str(i)]);
end