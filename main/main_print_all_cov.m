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


%% Setting
% hyperparameters
ell = 0.5;
sf = 1.5;
hyp = log([ell, sf]);

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
% x  = scale*rand(n-1, d);     x   = [x;   x(end, :)];                % duplication
% xd = scale*rand(nd-2, d);    xd  = [xd;  x(end, :);	xd(end, :)];	% duplication
% z  = scale*rand(ns-2, d);    z   = [z;   x(end, :);  xd(end, :)];  	% duplication

x = [
0.789963029944531, 0.111705744193203, 0.189710406017580;
0.318524245398992, 0.136292548938299, 0.495005824990221;
0.534064127370726, 0.678652304800188, 0.147608221976689;
0.089950678770581, 0.495177019089661, 0.054974146906188;
0.089950678770581, 0.495177019089661, 0.054974146906188];

xd = [
0.850712674289007, 0.929608866756663, 0.582790965175840;
0.560559527354885, 0.696667200555228, 0.815397211477421;
0.089950678770581, 0.495177019089661, 0.054974146906188;
0.560559527354885, 0.696667200555228, 0.815397211477421];

z = [
0.879013904597178, 0.988911616079589, 0.000522375356945;
0.089950678770581, 0.495177019089661, 0.054974146906188;
0.560559527354885, 0.696667200555228, 0.815397211477421];

xxd = [zeros(n, 1),   x;
       1*ones(nd, 1), xd;
       2*ones(nd, 1), xd;
       3*ones(nd, 1), xd];

% We have function and derivative observation at derivative inputs.
% Or we have surface normal vectors at sampled hit points.
% x = [x; xd];

print_matrix_for_reference(x, 'x');
print_matrix_for_reference(xd, 'xd');
print_matrix_for_reference(z, 'z');


%% Covariances
print_matrix_for_reference(sq_dist(x'),     'sqDist(x)');
print_matrix_for_reference(sq_dist(x', z'), 'sqDist(x, z)');
print_matrix_for_reference(delta(x, x, 1),  'delta(x1, x1)');
print_matrix_for_reference(delta(x, z, 1),  'delta(x1, z1)');

print_cov({@covSEiso}, 'covSEiso', hyp, x, z);
print_covDerObs(@covSEisoDerObs, 'covSEiso', hyp, xxd, z);

print_cov({@covMaterniso, 3}, 'covMaterniso3', hyp, x, z);
print_covDerObs(@covMaterniso3DerObs, 'covMaterniso3', hyp, xxd, z);

print_cov({@covSparseiso}, 'covSparseiso', hyp, x, z);
print_covDerObs(@covSparseisoDerObs, 'covSparseiso', hyp, xxd, z);

