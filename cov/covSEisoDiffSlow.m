function k = covSEisoDiffSlow(hyp, x, z, i, pdx, pdz)

%% function name convention
% cov:  covariance function
% SE:   squared exponential
% iso:  isotropic
% Diff: differentiable w.r.t input coordinates or take derivative observations
% Slow: non-vectorized or element-wised version

%% input/output arguments
% hyp:  [1x2]   hyperparameters, hyp = [log(ell), log(sigma_f)]
% x:    [dx1]   first function/derivative input vector
% z:    [dx1]   second function/derivative input vector, default: [] meaning z = x
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
   
% default parameters
if nargin < 2, K = '2'; return; end                  % report number of parameters
if nargin < 3, z = [];  end                                  % make sure, z exists
xeqz = numel(z)==0; dg = strcmp(z,'diag') && numel(z)>0;        % determine mode
if nargin < 4, i = 0;   end
if nargin < 5, pdx = 0; end
if nargin < 5, pdz = 0; end

% component function
f_handles.k = @k;
% no derivative w.r.t hyperparameters
if i ~= 0
    if pdx == 0
        if pdz == 0
            % k(x, z)
            k = 
        else
            % k(x, dz/dz_j)
        end
    else
        if pdz == 0
            % k(dx/dx_i, z)
        else
            % k(dx/dx_j, dz/dz_j)
        end
    end
    
% derivative w.r.t hyperparameters
else
    if i == 1
        if pdx == 0
            if pdz == 0
                % dk(x, z)/dlog(ell)
            else
                % dk(x, dz/dz_j)/dlog(ell)
            end
        else
            if pdz == 0
                % dk(dx/dx_i, z)/dlog(ell)
            else
                % dk(dx/dx_j, dz/dz_j)/dlog(ell)
            end
        end
    else
        if pdx == 0
            if pdz == 0
                % dk(x, z)/dlog(sigma_f)
            else
                % dk(x, dz/dz_j)/dlog(sigma_f)
            end
        else
            if pdz == 0
                % dk(dx/dx_i, z)/dlog(sigma_f)
            else
                % dk(dx/dx_j, dz/dz_j)/dlog(sigma_f)
            end
        end
    end
end

%% component functions
function dk_ds()
end

	*					which requires four components, ds/dx_i, ds/dx'_j, d^2k/ds^2 and d^2s/dx_i dx'_j
	*					as well as dk/ds which is already required in Isotropic.
	*

% Copyright (c) by Soohwan Kim, 2014-06-11.