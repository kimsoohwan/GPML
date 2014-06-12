clc
clear all
close all

addpath('../util/');
addpath('../cov/');
addpath('../cov/covSEisoDiff');

% hyperparameters
hyp = log([1, 1]);

% data
d = 3;
n  = 5;  x  = rand(n, d);
nd = 4;  xd = rand(nd, d);
ns = 3;  z  = rand(ns, d);

%% Test: covSEisoDiffCwise vs. covSEiso
testcase_name = '[covSEisoDiffCwise vs. covSEiso]';
disp(testcase_name);

% K(x, x)
TEST_EQ(covSEisoDiffCwise(hyp, x), ...
        covSEiso(hyp, x), ...
        '\tK(x, x)');

% dK(x, x)/dtheta_i
for i = 1:length(hyp)
    TEST_EQ(covSEisoDiffCwise(hyp, x, [], i), ...
            covSEiso(hyp, x, [], i), ...
            ['\tdK(x, x)/dtheta_', num2str(i)]);
end

% K(x, z)
TEST_EQ(covSEisoDiffCwise(hyp, x, z), ...
        covSEiso(hyp, x, z), ...
        '\tK(x, z)');


%% Test: covSEisoDiffBwise vs. covSEiso
testcase_name = '[covSEisoDiffBwise vs. covSEiso]';
disp(testcase_name);

% K(x, x)
TEST_EQ(covSEisoDiffBwise(hyp, x), ...
        covSEiso(hyp, x), ...
        '\tK(x, x)');

% dK(x, x)/dtheta_i
for i = 1:length(hyp)
    TEST_EQ(covSEisoDiffBwise(hyp, x, [], i), ...
            covSEiso(hyp, x, [], i), ...
            ['\tdK(x, x)/dtheta_', num2str(i)]);
end

% K(x, z)
TEST_EQ(covSEisoDiffBwise(hyp, x, z), ...
        covSEiso(hyp, x, z), ...
        '\tK(x, z)');


%% Test: covSEisoDiffCwise vs. covSEisoDiffBwise
testcase_name = '[covSEisoDiffCwise vs. covSEisoDiffBwise]';
disp(testcase_name);

% K(x, x)
TEST_EQ(covSEisoDiffCwise(hyp, x, [], 0, xd), ...
        covSEisoDiffBwise(hyp, x, [], 0, xd), ...
        '\tK(x, x)');

% dK(x, x)/dtheta_i    
for i = 0:length(hyp)                     
    TEST_EQ(covSEisoDiffCwise(hyp, x, [], i, xd), ...
            covSEisoDiffBwise(hyp, x, [], i, xd), ...
            ['\tdK(x, x)/dtheta_', num2str(i)]);
end

% K(x, z)
TEST_EQ(covSEisoDiffCwise(hyp, x, z, 0, xd), ...
        covSEisoDiffBwise(hyp, x, z, 0, xd), ...
        '\tK(x, z)');