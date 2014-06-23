function K = covMaterniso3Diff(hyp, x, z, i, xd, f_Bwise)

%% function name convention
% cov:      covariance function
% Matern:   Matern
% iso:      isotropic
% Diff:     differentiable w.r.t input coordinates or take derivative observations

%% input/output arguments
% hyp:  [1x2]   hyperparameters, hyp = [log(ell), log(sigma_f)]
% x:    [nxd]   first function input vectors
% z:    [nsxd]  second function input vectors, default: [] meaning z = x
% i:            partial deriavtive coordiante w.r.t hyperparameters, default: 0
% xd:   [ndxd]  first derivative input vectors
% f_Bwise:      true: block matrix-wised version, false: coeffient-wised version
% K:    [nnxns]  covariance, nn = n + nd*d

%% default parameters
if nargin < 2, K = '2'; return; end                  % report number of parameters
if nargin < 3, z = [];  end                                   % make sure, z exists
if nargin < 4, i = 0;   end
if nargin < 5, xd = []; end
if nargin < 6, f_Bwise = true; end

%% component functions

% k = sigma_f^2 * (1-s) * exp(s)
f_handles.k             = @(hyp_, x_, z_, i_, pdx_, pdz_) exp(2*hyp_(2))*(1-s(hyp_, x_, z_)).*exp(s(hyp_, x_, z_));

% dk/ds   = k -   sigma_f^2 * exp(s)
% d2k/ds2 = k - 2*sigma_f^2 * exp(s)
% d3k/ds3 = k - 3*sigma_f^2 * exp(s)
f_handles.dk_ds         = @(hyp_, x_, z_, i_, pdx_, pdz_) f_handles.k(hyp_, x_, z_, i_, pdx_, pdz_) -   exp(2*hyp_(2))*exp(s(hyp_, x_, z_));
f_handles.d2k_ds2       = @(hyp_, x_, z_, i_, pdx_, pdz_) f_handles.k(hyp_, x_, z_, i_, pdx_, pdz_) - 2*exp(2*hyp_(2))*exp(s(hyp_, x_, z_));
f_handles.d3k_ds3       = @(hyp_, x_, z_, i_, pdx_, pdz_) f_handles.k(hyp_, x_, z_, i_, pdx_, pdz_) - 3*exp(2*hyp_(2))*exp(s(hyp_, x_, z_));

% ds/dxi =  3*(xi - zi)/(ell^2 * s)
% ds/dzj = -3*(xj - zj)/(ell^2 * s)
% d2s/dxi dzj = -3*d(i, j)/(ell^2 * s) - (1/s) * ds/dxi * ds/dzj
f_handles.ds_dxi        = @(hyp_, x_, z_, i_, pdx_, pdz_)  3*delta(x_, z_, pdx_)./(exp(2*hyp_(1))*s(hyp_, x_, z_));
f_handles.ds_dzj        = @(hyp_, x_, z_, i_, pdx_, pdz_) -3*delta(x_, z_, pdz_)./(exp(2*hyp_(1))*s(hyp_, x_, z_));
f_handles.d2s_dxi_dzj   = @(hyp_, x_, z_, i_, pdx_, pdz_)      -3*(pdx_ == pdz_)./(exp(2*hyp_(1))*s(hyp_, x_, z_)) ...
                                                               - f_handles.ds_dxi(hyp_, x_, z_, i_, pdx_, pdz_) ...
                                                               .*f_handles.ds_dzj(hyp_, x_, z_, i_, pdx_, pdz_) ...
                                                               ./s(hyp_, x_, z_);

% ds/dell           = (-1/ell)*s
% d2s/dell dxi      = (-1/ell)*ds/dxi
% d2s/dell dzj      = (-1/ell)*ds/dzj
% d3s/dell dxi dzj	= (-1/ell)*d2s/dxi dzj
f_handles.ds_dell           = @(hyp_, x_, z_, i_, pdx_, pdz_) (-1/exp(hyp_(1)))*s(hyp_, x_, z_);
f_handles.d2s_dell_dxi      = @(hyp_, x_, z_, i_, pdx_, pdz_) (-1/exp(hyp_(1)))*f_handles.ds_dxi(hyp_, x_, z_, i_, pdx_, pdz_);
f_handles.d2s_dell_dzj      = @(hyp_, x_, z_, i_, pdx_, pdz_) (-1/exp(hyp_(1)))*f_handles.ds_dzj(hyp_, x_, z_, i_, pdx_, pdz_);
f_handles.d3s_dell_dxi_dzj  = @(hyp_, x_, z_, i_, pdx_, pdz_) (-1/exp(hyp_(1)))*f_handles.d2s_dxi_dzj(hyp_, x_, z_, i_, pdx_, pdz_);

% call
if f_Bwise
    K = covisoDiffBwise(f_handles, hyp, x, xd, z, i);
else
    K = covisoDiffCwise(f_handles, hyp, x, xd, z, i);
end

%% sub component function
% s = -sqrt(3*r^2/ell^2)
function value = s(hyp_, x_, z_)
    ell = exp(hyp_(1));                             % ell
    value = -sqrt(3*sq_dist(x_'/ell, z_'/ell));     % s = -sqrt(3*r^2/ell^2)
end

end