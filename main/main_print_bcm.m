clc
clear all
close all

%% Setting
% hyperparameters
ell = 0.5;
sf = 1.5;
hyp = log([ell, sf]);

%% initialization
Sigma = zeros(4, 4);
Var = zeros(4, 1);
mu_by_Sigma = zeros(4, 1);
mu_by_Var = zeros(4, 1);

%% prediction one
mu1 = [0.0344    0.4387    0.3816    0.7655]';
x1 = [0.9134    0.5469    0.9706    0.1419;
      0.6324    0.9575    0.9572    0.4218;
      0.0975    0.9649    0.4854    0.9157;
      0.2785    0.1576    0.8003    0.7922];
Sigma1 = feval(@covSEiso, hyp, x1);
Var1 = diag(Sigma1);      
print_matrix_for_reference(mu1, 'mu1');
print_matrix_for_reference(Sigma1, 'Sigma1');
print_matrix_for_reference(Var1, 'Var1');

% update
Sigma = Sigma + inv(Sigma1);
mu_by_Sigma = mu_by_Sigma + Sigma1\mu1;
print_matrix_for_reference(Sigma, 'sum_inv_covariances 1');
print_matrix_for_reference(mu_by_Sigma, 'sum_weighted_means_by_cov 1');

Var = Var + 1./Var1;
mu_by_Var = mu_by_Var + (1./Var1).*mu1;
print_matrix_for_reference(Var, 'sum_inv_variances 1');
print_matrix_for_reference(mu_by_Var, 'sum_weighted_means_by_var 1');

%% prediction two
mu2 = [0.7952    0.1869    0.4898    0.4456]';
x2 = [0.6463    0.6797    0.4984    0.2238;
      0.7094    0.6551    0.9597    0.7513;
      0.7547    0.1626    0.3404    0.2551;
      0.2760    0.1190    0.5853    0.5060];
Sigma2 = feval(@covSEiso, hyp, x2);  
Var2 = diag(Sigma2);      
print_matrix_for_reference(mu2, 'mu2');
print_matrix_for_reference(Sigma2, 'Sigma2');
print_matrix_for_reference(Var2, 'Var2');

% update
Sigma = Sigma + inv(Sigma2);
mu_by_Sigma = mu_by_Sigma + Sigma2\mu2;
print_matrix_for_reference(Sigma, 'sum_inv_covariances 2');
print_matrix_for_reference(mu_by_Sigma, 'sum_weighted_means_by_cov 2');

Var = Var + 1./Var2;
mu_by_Var = mu_by_Var + (1./Var2).*mu2;
print_matrix_for_reference(Var, 'sum_inv_variances 2');
print_matrix_for_reference(mu_by_Var, 'sum_weighted_means_by_var 2');

%% prediction three
mu3 = [0.4733    0.3517    0.8308    0.5853]';
x3 = [0.6991    0.1386    0.2543    0.3500;
      0.8909    0.1493    0.8143    0.1966;
      0.9593    0.2575    0.2435    0.2511;
      0.5472    0.8407    0.9293    0.6160];
Sigma3 = feval(@covSEiso, hyp, x3);  
Var3 = diag(Sigma3);      
print_matrix_for_reference(mu3, 'mu3');
print_matrix_for_reference(Sigma3, 'Sigma3');
print_matrix_for_reference(Var3, 'Var3');

% update
Sigma = Sigma + inv(Sigma3);
mu_by_Sigma = mu_by_Sigma + Sigma3\mu3;
print_matrix_for_reference(Sigma, 'sum_inv_covariances 3');
print_matrix_for_reference(mu_by_Sigma, 'sum_weighted_means_by_cov 3');

Var = Var + 1./Var3;
mu_by_Var = mu_by_Var + (1./Var3).*mu3;
print_matrix_for_reference(Var, 'sum_inv_variances 3');
print_matrix_for_reference(mu_by_Var, 'sum_weighted_means_by_var 3');

%% final
Sigma = inv(Sigma);
mu_by_Sigma = Sigma*mu_by_Sigma;
print_matrix_for_reference(Sigma, 'Sigma');
print_matrix_for_reference(mu_by_Sigma, 'mu_by_Sigma');

Var = 1./Var;
mu_by_Var = Var.*mu_by_Var;
print_matrix_for_reference(Var, 'Var');
print_matrix_for_reference(mu_by_Var, 'mu_by_Var');