function k = covSEisoDiffBase(hyp, x, z, i, pdx, pdz)

%% function name convention
% cov:  covariance function
% SE:   squared exponential
% iso:  isotropic
% Diff: differentiable w.r.t input coordinates or take derivative observations
% Base: base function for covIsoDiff

%% input/output arguments
% hyp:  [1x2]   hyperparameters, hyp = [log(ell), log(sigma_f)]
% x:    [1xd]   first function/derivative input vector
% z:    [1xd]   second function/derivative input vector, default: [] meaning z = x
% i:            partial deriavtive coordiante w.r.t hyperparameters, default: 0
% pdx:          partial derivative coordinate w.r.t the first input
%               if dxi = 0 (default), x = function input vector, else x = derivative input
% pdz:          partial derivative coordinate w.r.t the second input
%               if pdz = 0 (default), z = function input vector, else z = derivative input
% k:            covariance

%% Squared Exponential covariance function with isotropic distance measure.
%
% It extends the original isotropic squared exponential covariance function
% to take derivative inputs as well as function inputs.
%
% The isotropic squared exponential covariance function is:
%
% k(x, x') = sigma_f^2 * f(s), f(s) = exp(s), s(r) = -r^2/(2*ell^2), r = |x-x'|
%
% (1) Partial derivatives with respective to input coordinates
%
%           dx           dk(x, x')     dk      ds
%   (a) k(------, dx') = ---------- = ---- * ------
%          dx_i             dx_i       ds     dx_i
%          
%               dx'       dk(x, x')     dk      ds
%   (b) k(dx, -------) = ---------- = ---- * --------
%              dx'_j        dx'_j       ds     dx'_j
%
%           dx      dx'       dk(x, x')      d^2k    ds    ds      dk     d^2s
%   (c) k(------, -------) = ------------ = ------*-----*------ + ----*-----------
%          dx_i    dx'_j      dx_i dx'_j     ds^2   dx_i  dx'_j    ds   dx_i dx'_j
%
%
% (2) Partial derivatives with respective to input coordinates and hyperparameters for learning
%
%              d
%    (0) -------------- k(?x, ?x') = 2 * k(?x, ?x')
%        dlog(sigma_f)
%
%            d         dx             d^2 k(x, x')           d^2k    ds     ds      dk     d^2s
%    (a) ---------- k(------, dx') = ---------------- = ell*(------*------*------ + ----*----------)
%         dlog(ell)    dx_i           dlog(ell) dx_i          ds^2   dell   dx_i     ds   dell dx_i 
%
%             d            dx'        d^2 k(x, x')           d^2k    ds     ds       dk     d^2s
%    (b) ---------- k(x, ------) = ----------------- = ell*(------*------*------- + ----*-----------)
%         dlog(ell)       dx'_j     dlog(ell) dx'_j          ds^2   dell   dx'_j     ds   dell dx'_j 
%
%              d         dx      dx'            d^3k    ds     ds    ds      d^2k     d^2s      ds       d^2k    ds      d^2s        d^2k    ds      d^2s        dk        d^3s
%     (c) ---------- k(------, -------) = ell*(------*------*-----*------ + ------*----------*------- + ------*-----*------------ + ------*------*----------- + ----*-----------------)
%          dlog(ell)    dx_i    dx'_j           ds^3   dell   dx_i  dx'_j    ds^2   dell dx_i  dx'_j     ds^2   dx_i  dell dx'_j     ds^2   dell   dx_i dx'_j    ds   dell dx_i dx'_j
%
%                                               d^3k    ds     ds    ds      d^2k      d^2s      ds        ds      d^2s          ds      d^2s         dk        d^3s
%                                       = ell*(------*------*-----*------ + ------*(----------*------- + ------*------------ + ------*-----------) + ----*-----------------)
%                                               ds^3   dell   dx_i  dx'_j    ds^2    dell dx_i  dx'_j     dx_i  dell dx'_j      dell   dx_i dx'_j     ds   dell dx_i dx'_j
   
% % default parameters
% % if nargin < 2, K = '2'; return; end                  % report number of parameters
% if nargin < 3, z = x;  end                                  % make sure, z exists
% % xeqz = numel(z)==0; dg = strcmp(z,'diag') && numel(z)>0;        % determine mode
% if nargin < 4, i = 0;   end
% if nargin < 5, pdx = 0; end
% if nargin < 6, pdz = 0; end

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
k = covisoDiff(f_handles, hyp, x, z, i, pdx, pdz);

%% sub component function
% s = -r^2/(2*ell^2)
function value = s(hyp_, x_, z_)
    ell = exp(hyp_(1));                             % ell
    value = (-1/2) * sq_dist(x_'/ell, z_'/ell);     % s = (-1/2)*r^2/ell^2 
end

end