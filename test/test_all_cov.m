clc
clear all
close all

%% path
addpath('../util/');

addpath('../cov/covSparseiso');

addpath('../cov/coviso');
addpath('../cov/coviso/covSEisoDiff');
addpath('../cov/coviso/covMaterniso3Diff');
addpath('../cov/coviso/covSparseisoDiff');

addpath('../cov/covDiff');
addpath('../cov/covDiff/covSEisoDiffFast');
addpath('../cov/covDiff/covMaterniso3DiffFast');
addpath('../cov/covDiff/covSparseisoDiffFast');

%% setting
% hyperparameters
hyp = log([1, 1]);

% data
d = 3;
% n  = 5;  x  = rand(n, d);
% nd = 4;  xd = rand(nd, d);
% ns = 3;  z  = rand(ns, d);
n  = 5;  x  = rand(n-1, d);     x   = [x;   x(end, :)];                 % duplication
nd = 4;  xd = rand(nd-2, d);    xd  = [xd;  x(end, :);	xd(end, :)];	% duplication
ns = 3;  z  = rand(ns-2, d);    z   = [z;   x(end, :);  xd(end, :)];  	% duplication


%% test cases
% testcase_covSEisoDiff;
% testcase_covMaterniso3Diff;
% testcase_covSparseisoDiff;

testcase_covSEisoDiffFast;
testcase_covMaterniso3DiffFast;
testcase_covSparseisoDiffFast;
