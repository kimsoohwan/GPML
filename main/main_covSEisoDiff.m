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

% We have function and derivative observation at derivative inputs.
% Or we have surface normal vectors at sampled hit points.
x = [x; xd];

% K
K = covSEisoDiffCwise(hyp, x, [], 0, xd);
print_matrix_for_reference(K);

% dK/dlog(theta)
dK_dlogell  = covSEisoDiffCwise(hyp, x, [], 1, xd);
dK_dlogsf   = covSEisoDiffCwise(hyp, x, [], 2, xd);
print_matrix_for_reference(dK_dlogell);
print_matrix_for_reference(dK_dlogsf);

% Ks
Ks = covSEisoDiffCwise(hyp, x, z, 0, xd);
print_matrix_for_reference(Ks);
