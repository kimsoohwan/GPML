function A = meanGP(hyp, x, i)
% GP mean function. The mean function does not have any parameters.

% global variables
global global_mean_func;
global global_cov_func;
global global_lik_func;
global global_inf_method;
global global_hyp;
global global_x;
global global_y;

if nargin<2, A = '0'; return; end             % report number of hyperparameters 

% % for now, it only works for 1D cases
% if size(x, 2) > 2
%     error('For now, it only works for 1D cases')
% end
% 
% if size(x, 2) == 2
%     mask = x(:, 1) == 0;
%     x = x(mask, 2);
% end
    
% mean
[dummy, dummy, A, dummy] = gp(global_hyp, global_inf_method, global_mean_func, global_cov_func, global_lik_func, global_x, global_y, x);
