clc
clear all
close all

%% path
addpath('../util/');

addpath('../cov/covSparseiso');

addpath('../cov/covisoDerObsUnstable');
addpath('../cov/covisoDerObsUnstable/covSEisoDerObsUnstable');
addpath('../cov/covisoDerObsUnstable/covMaterniso3DerObsUnstable');
addpath('../cov/covisoDerObsUnstable/covSparseisoDerObsUnstable');

addpath('../cov/covDerObs');
addpath('../cov/covDerObs/covSEisoDerObs');
addpath('../cov/covDerObs/covMaterniso3DerObs');
addpath('../cov/covDerObs/covSparseisoDerObs');

%% setting
% hyperparameters
hyp = log([0.2, 0.5]);

% data
d = 3;
n  = 5;
nd = 4;
ns = 3;
% scale = 0.1;
scale = 1;
% scale = 10;

% x  = scale*rand(n, d);
% xd = scale*rand(nd, d);
% z  = scale*rand(ns, d);
%
x	 = [0.0423    0.6751    0.0839
        0.9730    0.3610    0.9748
        0.1892    0.6203    0.6513
        0.6671    0.8112    0.2312
        0.5864    0.0193    0.4035];
xd   = [0.1220    0.1522    0.0943
        0.2684    0.3480    0.9300
        0.2578    0.1217    0.3990
        0.3317    0.8842    0.0474];
z    = [0.3424    0.5449    0.0548
        0.7360    0.6862    0.3037
        0.7947    0.8936    0.0462];  
       
% x  = scale*rand(n-1, d);     x   = [x;   x(end, :)];                 % duplication
% xd = scale*rand(nd-2, d);    xd  = [xd;  x(end, :);	xd(end, :)];	% duplication
% z  = scale*rand(ns-2, d);    z   = [z;   x(end, :);  xd(end, :)];  	% duplication

xx  = [zeros(n, 1), x];
xxd = [zeros(n, 1), x
       ones(nd, 1), xd];
    
% sigma_n
% sigma_n = 0;
sigma_n = 0.0000001;


%% test cases
% testcase_cov({@covSEiso},        'covSEiso',      hyp, x, sigma_n);
% testcase_cov({@covMaterniso, 3}, 'covMaterniso3', hyp, x, sigma_n);
% testcase_cov({@covSparseiso},    'covSparseiso',  hyp, x, sigma_n);

% testcase_covisoDerObsUnstable(@covSEisoDerObsUnstable,      {@covSEiso},        'covSEiso',      hyp, x, z, xx, xxd, sigma_n);
% testcase_covisoDerObsUnstable(@covMaterniso3DerObsUnstable, {@covMaterniso, 3}, 'covMaterniso3', hyp, x, z, xx, xxd, sigma_n);
% testcase_covisoDerObsUnstable(@covSparseisoDerObsUnstable,  {@covSparseiso},    'covSparseiso',  hyp, x, z, xx, xxd, sigma_n);

testcase_covDerObs(@covSEisoDerObs,      @covSEisoDerObsUnstable,      {@covSEiso},        'covSEiso',      hyp, x, z, xx, xxd, sigma_n);
testcase_covDerObs(@covMaterniso3DerObs, @covMaterniso3DerObsUnstable, {@covMaterniso, 3}, 'covMaterniso3', hyp, x, z, xx, xxd, sigma_n);
testcase_covDerObs(@covSparseisoDerObs,  @covSparseisoDerObsUnstable,  {@covSparseiso},    'covSparseiso',  hyp, x, z, xx, xxd, sigma_n);
