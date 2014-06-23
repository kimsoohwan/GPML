clc
clear all
close all

%% path
addpath('../util/');
addpath('../cov/');

addpath('../cov/coviso');
addpath('../cov/coviso/covSEisoDiff');
addpath('../cov/coviso/covMaterniso3Diff');
% addpath('../cov/coviso/covSpaseisoDiff');

addpath('../cov/covDiff');
addpath('../cov/covDiff/covSEisoDiffFast');
addpath('../cov/covDiff/covMaterniso3DiffFast');

%% test cases
% testcase_covSEisoDiff;
% testcase_covMaterniso3Diff;
% testcase_covSparseisoDiff;

% testcase_covSEisoDiffFast;
testcase_covMaterniso3DiffFast;

% TODO
% 1. covSparseiso
% 2. testcase_covSparseisoDiff;
% 3. covDerivative¶û ºñ±³ÇÏ±â
