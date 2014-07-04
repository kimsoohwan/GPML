function K = covSEisoDerObsUnstable(hyp, x, z, i, f_Bwise)

%% function name convention
% cov:      covariance function
% SE:       squared exponential
% iso:      isotropic
% DerObs:     differentiable w.r.t input coordinates or take derivative observations
% Unstable: cannot handle divide-by-zero

%% input/output arguments
% hyp:  [1x2]       hyperparameters, hyp = [log(ell), log(sigma_f)]
% x:    [nx(d+1)]   first function/derivative input vectors
% z:    [nsxd]      second function input vectors, default: [] meaning z = x
% i:                partial deriavtive coordiante w.r.t hyperparameters, default: 0
% f_Bwise:          true: block matrix-wised version, false: coeffient-wised version
% K:    [nnxns]     covariance, nn = n + nd*d

%% default parameters
if nargin < 2, K = '2'; return; end                  % report number of parameters
if nargin < 3, z = [];  end                                   % make sure, z exists
if nargin < 4, i = 0;   end
if nargin < 5, f_Bwise = true; end
dg = strcmp(z,'diag') && numel(z)>0;    % determine mode

% x, xd
if dg
    xd = [];
else
    % assume derivative observations are repeated d times
    d = size(x, 2) - 1; % first column = index
    xd_old = [];
    for dd = 1:d
        mask = x(:, 1) == d;
        xd = x(mask, 2:end);
        if dd > 1            
            assert(~any(xd_old(:) - xd(:)));
        end
        xd_old = xd;
    end
    
    xd = x(x(:, 1) == 1, 2:end);
    x  = x(x(:, 1) == 0, 2:end);
end

%% component functions

% k = sigma_f^2 * exp(s)
f_handles.k             = @(hyp_, x_, z_, i_, pdx_, pdz_) exp(2*hyp_(2))*exp(s(hyp_, x_, z_));

% dk/ds = d2k/ds2 = d3k/ds3 = k
f_handles.dk_ds         = f_handles.k;
f_handles.d2k_ds2       = f_handles.k;
f_handles.d3k_ds3       = f_handles.k;

% ds/dxi =   -(xi - zi)/ell^2
% ds/dzj =    (xj - zj)/ell^2
% d2s/dxi dzj = d(i, j)/ell^2
f_handles.ds_dxi        = @(hyp_, x_, z_, i_, pdx_, pdz_) -delta(x_, z_, pdx_)/exp(2*hyp_(1));
f_handles.ds_dzj        = @(hyp_, x_, z_, i_, pdx_, pdz_)  delta(x_, z_, pdz_)/exp(2*hyp_(1));
f_handles.d2s_dxi_dzj   = @(hyp_, x_, z_, i_, pdx_, pdz_)       (pdx_ == pdz_)/exp(2*hyp_(1));

% ds/dell           = (-2/ell)*s
% d2s/dell dxi      = (-2/ell)*ds/dxi
% d2s/dell dzj      = (-2/ell)*ds/dzj
% d3s/dell dxi dzj	= (-2/ell)*d2s/dxi dzj
f_handles.ds_dell           = @(hyp_, x_, z_, i_, pdx_, pdz_) (-2/exp(hyp_(1)))*s(hyp_, x_, z_);
f_handles.d2s_dell_dxi      = @(hyp_, x_, z_, i_, pdx_, pdz_) (-2/exp(hyp_(1)))*f_handles.ds_dxi(hyp_, x_, z_, i_, pdx_, pdz_);
f_handles.d2s_dell_dzj      = @(hyp_, x_, z_, i_, pdx_, pdz_) (-2/exp(hyp_(1)))*f_handles.ds_dzj(hyp_, x_, z_, i_, pdx_, pdz_);
f_handles.d3s_dell_dxi_dzj  = @(hyp_, x_, z_, i_, pdx_, pdz_) (-2/exp(hyp_(1)))*f_handles.d2s_dxi_dzj(hyp_, x_, z_, i_, pdx_, pdz_);

% call
if f_Bwise
    K = covisoDerObsBwiseUnstable(f_handles, hyp, x, xd, z, i);
else
    K = covisoDerObsCwiseUnstable(f_handles, hyp, x, xd, z, i);
end

%% sub component function
% s = -r^2/(2*ell^2)
function value = s(hyp_, x_, z_)
    ell = exp(hyp_(1));                             % ell
    value = (-1/2) * sq_dist(x_'/ell, z_'/ell);     % s = (-1/2)*r^2/ell^2 
end

end