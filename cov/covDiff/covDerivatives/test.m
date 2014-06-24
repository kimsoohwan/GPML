clc
clear all
close all

format long

%% K
% X = [2     2     2;
%      8     9     6;
%      2     3     4;
%      8     1     3];
%   
% Xs = [8     3     1;
%       5     5     5;
%       5     0     4;
%       9     0     0;
%       2     5     3;
%       7     7     1;
%       7     9     7];

X = [0.118997681558377,   0.223811939491137,   0.890903252535799;
	 0.498364051982143,   0.751267059305653,   0.959291425205444;
     0.959743958516081,   0.255095115459269,   0.547215529963803;
     0.340385726666133,   0.505957051665142,   0.138624442828679;
     0.585267750979777,   0.699076722656686,   0.149294005559057];
 
Xs = [0.257508254123736,   0.243524968724989,   0.251083857976031;
	  0.840717255983663,   0.929263623187228,   0.616044676146639;
      0.254282178971531,   0.349983765984809,   0.473288848902729;
      0.814284826068816,   0.196595250431208,   0.351659507062997];                    
  
% squared distance
sqDist11 = sq_dist(X');
sqDist12 = sq_dist(X', Xs');

% setting
inf_method  = @infExact;
mean_func   = {@meanZero};
lik_func    = {@likGauss};

hyp.mean = [];
% hyp.cov = log([1.5, 2.5]);
% hyp.lik = log(0.3);
hyp.cov = log([0.5, 1.5]);
hyp.lik = log(0.1);
covFunc = {@covSEiso};
% covFunc = {@covMaterniso, 3};

K       = feval(covFunc{:}, hyp.cov, X);
K_ell   = feval(covFunc{:}, hyp.cov, X, [], 1);
K_sf2   = feval(covFunc{:}, hyp.cov, X, [], 2);
Ks     	= feval(covFunc{:}, hyp.cov, X, Xs);
Kss     = feval(covFunc{:}, hyp.cov, Xs, 'diag');

% inference
% y = [1;
%      4;
%      3;
%      7];
y = [0.913337361501670;
     0.152378018969223;
     0.825816977489547;
     0.538342435260057;
     0.996134716626885];
 
[y_mu, y_s2, f_mu, f_s2] = gp(hyp, inf_method, mean_func, covFunc, lik_func, X, y, Xs);

% partial derivatives of nlZ
[post, nlZ, dnlZ] = feval(inf_method, hyp, mean_func, covFunc, lik_func, X, y);

% training
hyp0 = hyp;
[hyp1, nlZ, iter] = minimize(hyp0, @gp, -50, inf_method, mean_func, covFunc, lik_func, X, y);

%% K_FDI

% covarnace between function and derivatives   
covFunc = {@covSEisoFD};
K_FD1       = feval(covFunc{:}, hyp.cov, X, X, 1);
K_FD2       = feval(covFunc{:}, hyp.cov, X, X, 2);
K_FD3       = feval(covFunc{:}, hyp.cov, X, X, 3);
K_FD1_ell   = feval(covFunc{:}, hyp.cov, X, X, 1, 1);
K_FD1_sf    = feval(covFunc{:}, hyp.cov, X, X, 1, 2);
K_FD2_ell   = feval(covFunc{:}, hyp.cov, X, X, 2, 1);
K_FD2_sf    = feval(covFunc{:}, hyp.cov, X, X, 2, 2);
K_FD3_ell   = feval(covFunc{:}, hyp.cov, X, X, 3, 1);
K_FD3_sf    = feval(covFunc{:}, hyp.cov, X, X, 3, 2);

covFunc = {@covSEisoDD};
K_DD11       = feval(covFunc{:}, hyp.cov, X, 1, X, 1);
K_DD12       = feval(covFunc{:}, hyp.cov, X, 1, X, 2);
K_DD13       = feval(covFunc{:}, hyp.cov, X, 1, X, 3);
K_DD21       = feval(covFunc{:}, hyp.cov, X, 2, X, 1);
K_DD22       = feval(covFunc{:}, hyp.cov, X, 2, X, 2);
K_DD23       = feval(covFunc{:}, hyp.cov, X, 2, X, 3);
K_DD31       = feval(covFunc{:}, hyp.cov, X, 3, X, 1);
K_DD32       = feval(covFunc{:}, hyp.cov, X, 3, X, 2);
K_DD33       = feval(covFunc{:}, hyp.cov, X, 3, X, 3);
K_DD11_ell       = feval(covFunc{:}, hyp.cov, X, 1, X, 1, 1);
K_DD12_ell       = feval(covFunc{:}, hyp.cov, X, 1, X, 2, 1);
K_DD13_ell       = feval(covFunc{:}, hyp.cov, X, 1, X, 3, 1);
K_DD21_ell       = feval(covFunc{:}, hyp.cov, X, 2, X, 1, 1);
K_DD22_ell       = feval(covFunc{:}, hyp.cov, X, 2, X, 2, 1);
K_DD23_ell       = feval(covFunc{:}, hyp.cov, X, 2, X, 3, 1);
K_DD31_ell       = feval(covFunc{:}, hyp.cov, X, 3, X, 1, 1);
K_DD32_ell       = feval(covFunc{:}, hyp.cov, X, 3, X, 2, 1);
K_DD33_ell       = feval(covFunc{:}, hyp.cov, X, 3, X, 3, 1);
K_DD11_sf        = feval(covFunc{:}, hyp.cov, X, 1, X, 1, 2);
K_DD12_sf        = feval(covFunc{:}, hyp.cov, X, 1, X, 2, 2);
K_DD13_sf        = feval(covFunc{:}, hyp.cov, X, 1, X, 3, 2);
K_DD21_sf        = feval(covFunc{:}, hyp.cov, X, 2, X, 1, 2);
K_DD22_sf        = feval(covFunc{:}, hyp.cov, X, 2, X, 2, 2);
K_DD23_sf        = feval(covFunc{:}, hyp.cov, X, 2, X, 3, 2);
K_DD31_sf        = feval(covFunc{:}, hyp.cov, X, 3, X, 1, 2);
K_DD32_sf        = feval(covFunc{:}, hyp.cov, X, 3, X, 2, 2);
K_DD33_sf        = feval(covFunc{:}, hyp.cov, X, 3, X, 3, 2);

n2 = 2;
% covFunc = {@covSEisoFDI, n2};
% covFunc = {@covMaterniso3FDI, n2};
covFunc = {@covSparseisoFDI, n2};
hyp.cov = log([20, 2.5]);

K_FDI       = feval(covFunc{:}, hyp.cov, X);
K_FDI_ell   = feval(covFunc{:}, hyp.cov, X, [], 1);
K_FDI_sf    = feval(covFunc{:}, hyp.cov, X, [], 2);
Ks_FDI      = feval(covFunc{:}, hyp.cov, X, Xs);

% covarnace between function and derivatives
X = [0     2     2     2;
     0     8     9     6;
     1     2     2     2;
     1     8     9     6;
     2     2     2     2;
     2     8     9     6;
     3     2     2     2;
     3     8     9     6;
     0     2     3     4;
     0     8     1     3];
   
% covFunc = {@covSEisoDerivatives};
% covFunc = {@covMaternisoDerivatives};
covFunc = {@covSparseisoDerivatives};
K_Derivatives       = feval(covFunc{:}, hyp.cov, X);
K_Derivatives_ell   = feval(covFunc{:}, hyp.cov, X, [], 1);
K_Derivatives_sf    = feval(covFunc{:}, hyp.cov, X, [], 2);
Ks_Derivatives      = feval(covFunc{:}, hyp.cov, X, Xs);

if any(abs(K_FDI(:) - K_Derivatives(:)) > 1e-7)
    error('K_FDI != K_Derivatives');
end
if any(abs(K_FDI_ell(:) - K_Derivatives_ell(:)) > 1e-7)
    warning('K_FDI_ell != K_Derivatives_ell');
end
if any(abs(K_FDI_sf(:) - K_Derivatives_sf(:)) > 1e-7)
    warning('K_FDI_sf != K_Derivatives_sf');
end
if any(abs(Ks_FDI(:) - Ks_Derivatives(:)) > 1e-7)
    warning('Ks_FDI != Ks_Derivatives');
end

% isnan
if any(isnan(K_Derivatives))
    warning('K_FDI is not PD!');
end
if any(isnan(K_Derivatives_ell))
    warning('K_FDI_ell is not PD!');
end
if any(isnan(K_Derivatives_sf))
    warning('K_FDI_sf is not PD!');
end
if any(isnan(Ks_Derivatives))
    warning('Ks_FDI is not PD!');
end

% positive definite
[R, p] = chol(K_Derivatives);
if p ~= 0
    warning('K_FDI is not PD!');
end
[R, p] = chol(K_Derivatives_ell);
if p ~= 0
    warning('K_FDI_ell is not PD!');
end
[R, p] = chol(K_Derivatives_sf);
if p ~= 0
    warning('K_FDI_sf is not PD!');
end

% setting
inf_method  = @infExactDerivatives;
mean_func   = {@meanZero};
lik_func    = {@likGaussDerivatives};

hyp.mean = [];
hyp.cov = log([1.5, 2.5]);
hyp.lik = log([0.3, 0.5]); % snf, snd
covFunc = {@covSEisoDerivatives};
covFunc = {@covMaternisoDerivatives};

y = [1;
     4;
     1;
     -1;
     2;
     -2;
     3;
     -3;
     3;
     7];

[y_mu, y_s2, f_mu, f_s2] = gp(hyp, inf_method, mean_func, covFunc, lik_func, X, y, Xs);

% partial derivatives of nlZ
[post, nlZ, dnlZ] = feval(inf_method, hyp, mean_func, covFunc, lik_func, X, y);

% training

hyp0 = hyp;
[hyp1, nlZ, iter] = minimize(hyp0, @gp, -50, inf_method, mean_func, covFunc, lik_func, X, y);
