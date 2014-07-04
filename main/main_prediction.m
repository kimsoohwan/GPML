clc
clear all
close all

%% path
addpath('../util/');

% mean function
addpath('../mean');

% covariance function
addpath('../cov/covSparseiso');

addpath('../cov/covisoDerObsUnstable');
addpath('../cov/covisoDerObsUnstable/covSEisoDerObsUnstable');
addpath('../cov/covisoDerObsUnstable/covMaterniso3DerObsUnstable');
addpath('../cov/covisoDerObsUnstable/covSparseisoDerObsUnstable');

addpath('../cov/covDerObs');
addpath('../cov/covDerObs/covSEisoDerObs');
addpath('../cov/covDerObs/covMaterniso3DerObs');
addpath('../cov/covDerObs/covSparseisoDerObs');

% likelihood function
addpath('../lik');

% inference method
addpath('../inf');


%% data
% function
f  = @(x) 2*sin(3*pi*x) - 2*x.^3 + x.^2 -3*x + 1 - sqrt(x) - cos(10*pi*x);
df = @(x) 6*pi*cos(3*pi*x) - 6*x.^2 + 2*x - 3 - 1./(2*sqrt(x)) + 10*pi*sin(10*pi*x);

% points
x = (0:0.01:1)';
y = feval(f, x);

% observations
xx = (0.04:0.1:0.95)';
% sn = 0.1; snd = 0.2;
% yy = feval(f, xx) + sn*randn(size(xx));
% dy = feval(df, xx) + snd*randn(size(xx));
yy = [1.162470824949019; 2.649518079414324; 0.826178856519631; -0.296449902131397; -2.925889279464160; -3.059479729867036; -2.684443970935226; -0.724953106178071; -0.871446982536995; -1.909489684161842];
dy = [41.704618124666638; -28.757549228476215; 14.122022625746855; -52.574386176465111; 15.885794730520876; -27.330344334593999; 43.308255238127209; -20.443387861589134; 22.877410410575429; -52.447400235726114];
 
% plot
figure('position', [0, 0, 1200, 400]);
subplot(1, 3, 1);
hold on
delta_x = 0.03;
plot(x, y, 'k-');
plot(xx, yy, 'or');
plot([xx' - delta_x; xx' + delta_x], [yy' - dy'*delta_x; yy' + dy'*delta_x], 'r-');
axis([0, 1, min(y)-0.1, max(y)+0.1]);

%% GP - Function Observations
% setting
mean_func = @meanZero;
cov_func = {@covSEiso};
% cov_func = {@covMaterniso, 3};
lik_func = @likGauss;
inf_method = @infExact;

% hyperparameter
ell = 0.5; sf = 1.5; sn = 0.1; 
hyp.cov = log([ell, sf]);
hyp.lik = log(sn);

% training
hyp1 = minimize(hyp, @gp, -100, inf_method, mean_func, cov_func, lik_func, xx, yy);

% prediction - regression
[dummy, dummy, fmu, fs2] = gp(hyp1, inf_method, mean_func, cov_func, lik_func, xx, yy, x);

% plot
subplot(1, 3, 2);
hold on;
ff = [fmu+2*sqrt(fs2); flipdim(fmu-2*sqrt(fs2),1)]; 
fill([x; flipdim(x,1)], ff, [7 7 7]/8)
plot(x, fmu, 'b-');
plot(x, y, 'k-');
plot(xx, yy, 'or');
axis([0, 1, min(y)-0.1, max(y)+0.1]);


%% GP - Derivative and Function Observations
% setting
mean_func = @meanZeroDerObs;
cov_func  = {@covSEisoDerObs};
lik_func  = @likGaussDerObs;
inf_method = @infExactDerObs;

% data
xxd = [zeros(length(xx), 1), xx;
       ones(length(xx), 1),  xx];
yyd = [yy;
       dy];
% xxd = [zeros(length(xx), 1), xx];
% yyd = [yy];

% hyperparameter
hyp.cov = log([ell, sf]);
hyp.lik = log([sn, sn]);

% training
hyp2 = minimize(hyp, @gp, -100, inf_method, mean_func, cov_func, lik_func, xxd, yyd);

% prediction - regression
[dummy, dummy, fmu, fs2] = gp(hyp2, inf_method, mean_func, cov_func, lik_func, xxd, yyd, x);

% plot
subplot(1, 3, 3);
hold on;
ff = [fmu+2*sqrt(fs2); flipdim(fmu-2*sqrt(fs2),1)]; 
fill([x; flipdim(x,1)], ff, [7 7 7]/8)
plot(x, fmu, 'b-');
plot(x, y, 'k-');
plot(xx, yy, 'or');
axis([0, 1, min(y)-0.1, max(y)+0.1]);