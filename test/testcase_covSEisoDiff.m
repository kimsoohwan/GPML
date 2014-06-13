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

% K(x, x)
TEST_EQ(covSEisoDiffCwise(hyp, x), ...
        covSEiso(hyp, x), ...
        'K(x, x)');

% dK(x, x)/dtheta_i
for i = 1:length(hyp)
    TEST_EQ(covSEisoDiffCwise(hyp, x, [], i), ...
            covSEiso(hyp, x, [], i), ...
            ['dK(x, x)/dtheta_', num2str(i)]);
end

% K(x, z)
TEST_EQ(covSEisoDiffCwise(hyp, x, z), ...
        covSEiso(hyp, x, z), ...
        'K(x, z)');


%% Test: covSEisoDiffBwise vs. covSEiso
test_name = '\t[TESE: covSEisoDiffBwise vs. covSEiso]\n';
fprintf(1, test_name);

% K(x, x)
TEST_EQ(covSEisoDiffBwise(hyp, x), ...
        covSEiso(hyp, x), ...
        'K(x, x)');

% dK(x, x)/dtheta_i
for i = 1:length(hyp)
    TEST_EQ(covSEisoDiffBwise(hyp, x, [], i), ...
            covSEiso(hyp, x, [], i), ...
            ['dK(x, x)/dtheta_', num2str(i)]);
end

% K(x, z)
TEST_EQ(covSEisoDiffBwise(hyp, x, z), ...
        covSEiso(hyp, x, z), ...
        'K(x, z)');


%% Test: covSEisoDiffCwise vs. covSEisoDiffBwise
test_name = '\t[TESE: covSEisoDiffCwise vs. covSEisoDiffBwise]\n';
fprintf(1, test_name);

% K(x, x)
TEST_EQ(covSEisoDiffCwise(hyp, x, [], 0, xd), ...
        covSEisoDiffBwise(hyp, x, [], 0, xd), ...
        'K(x, x)');

% dK(x, x)/dtheta_i    
for i = 1:length(hyp)                     
    TEST_EQ(covSEisoDiffCwise(hyp, x, [], i, xd), ...
            covSEisoDiffBwise(hyp, x, [], i, xd), ...
            ['dK(x, x)/dtheta_', num2str(i)]);
end

% K(x, z)
TEST_EQ(covSEisoDiffCwise(hyp, x, z, 0, xd), ...
        covSEisoDiffBwise(hyp, x, z, 0, xd), ...
        'K(x, z)');
    
%% Test: covSEisoDiffBwise - symmetric
test_name = '\t[TESE: covSEisoDiffBwise - symmetric]\n';
fprintf(1, test_name);

% K(x, x)
K = covSEisoDiffBwise(hyp, x, [], 0, xd);
TEST_EQ(K, ...
        K', ...
        'K(x, x)');

% dK(x, x)/dtheta_i    
for i = 1:length(hyp)
    K = covSEisoDiffBwise(hyp, x, [], i, xd);
    TEST_EQ(K, ...
            K', ...
            ['dK(x, x)/dtheta_', num2str(i)]);
end

%% Test: covSEisoDiffBwise - positive definite
test_name = '\t[TESE: covSEisoDiffBwise - positive definite]\n';
fprintf(1, test_name);

% sigma_n
sigma_n = 1;

% K(x, x)
K = covSEisoDiffBwise(hyp, x, [], 0, xd);
K = K + sigma_n^2*eye(size(K));
[R, p] = chol(K);
TEST_EQ(p, ...
        0', ...
        'K(x, x)');

% dK(x, x)/dtheta_i    
for i = 1:length(hyp)
    K = covSEisoDiffBwise(hyp, x, [], i, xd);
    K = K + sigma_n^2*eye(size(K));
    [R, p] = chol(K);
    TEST_EQ(p, ...
            0', ...
            ['dK(x, x)/dtheta_', num2str(i)]);
end