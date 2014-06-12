function k = covIsoDiff(hyp, x, z, i, pdx, pdz, f_handles)

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
% f_handles:    component function handles

%% Squared Exponential covariance function with isotropic distance measure.
%
% It extends the original isotropic squared exponential covariance function
% to take derivative inputs as well as function inputs.
%
% The isotropic squared exponential covariance function is:
%
% k(x, x') = sigma_f^2 * f(s), f(s) = exp(s), s(r) = -r^2/(2*ell^2), r = |x-x'|
%
% (1) Partial derivatives with respective to hyperparameters
%          
%             dk              dk          dsigma_f
%   (0) --------------- = ---------- * -------------- = 2 * k
%        dlog(sigma_f)     dsigma_f    dlog(sigma_f)
%
%           dk          dk       dell              dk      ds
%   (a) ----------- = ----- * ----------- = ell * ---- * ------
%        dlog(ell)     dell    dlog(ell)           ds     dell
%
%
% (2) Partial derivatives with respective to input coordinates
%
%           dx          dk(x, z)     dk      ds
%   (a) k(------, z) = ---------- = ---- * ------
%          dx_i           dx_i       ds     dx_i
%          
%              z        dk(x, z)     dk      ds
%   (b) k(x, ------) = ---------- = ---- * ------
%             dz_j        dz_j       ds     dz_j
%
%           dx       dz       dk(x, z)      d^2k    ds    ds      dk     d^2s
%   (c) k(------, -------) = ----------- = ------*-----*------ + ----*-----------
%          dx_i     dz_j      dx_i dz_j     ds^2   dx_i  dz_j     ds   dx_i dz_j
%
%
% (3) Partial derivatives with respective to input coordinates and hyperparameters for learning
%
%              d
%    (0) -------------- k(?x, ?z) = 2 * k(?x, ?z)
%        dlog(sigma_f)
%
%            d          dx            d^2 k(x, z)           d^2k    ds     ds      dk     d^2s
%    (a) ---------- k(------, z) = ---------------- = ell*(------*------*------ + ----*-----------)
%         dlog(ell)    dx_i         dlog(ell) dx_i          ds^2   dell   dx_i     ds   dell dx_i 
%
%             d            dz        d^2 k(x, z)            d^2k    ds     ds      dk     d^2s
%    (b) ---------- k(x, ------) = ---------------- = ell*(------*------*------ + ----*-----------)
%         dlog(ell)       dz_j      dlog(ell) dz_j          ds^2   dell   dz_j     ds   dell dz_j 
%
%             d         dx      dz            d^3k    ds     ds    ds      d^2k      d^2s      ds      d^2k    ds      d^2s        d^2k    ds      d^2s        dk        d^3s
%    (c) ---------- k(------, ------) = ell*(------*------*-----*------ + ------*-----------*------ + ------*------*----------- + ------*------*----------- + ----*----------------)
%         dlog(ell)    dx_i    dz_j           ds^3   dell   dx_i  dz_j     ds^2   dell dx_i   dz_j     ds^2   dx_i   dell dz_j     ds^2   dell   dx_i dz_j     ds   dell dx_i dz_j
%
%                                              d^3k    ds     ds    ds      d^2k      d^2s       ds       ds      d^2s         ds      d^2s         dk        d^3s
%                                      = ell*(------*------*-----*------ + ------*(-----------*------ + ------*----------- + ------*-----------) + ----*----------------)
%                                              ds^3   dell   dx_i  dz_j     ds^2    dell dx_i   dz_j     dx_i   dell dz_j     dell   dx_i dz_j      ds   dell dx_i dz_j

% hyperparameters
ell = exp(hyp(1));                                 % characteristic length scale

% no derivative w.r.t hyperparameters
if i == 0
    if pdx == 0
        if pdz == 0
            % [0]  k(x, z)
            k = f_handles.k(hyp, x, z, i, pdx, pdz);
        else
            % [2.b] k(x, dz/dz_j)
            k = f_handles.dk_ds(hyp, x, z, i, pdx, pdz) ...
                .* f_handles.ds_dzj(hyp, x, z, i, pdx, pdz);
        end
    else
        if pdz == 0
            % [2.a] k(dx/dx_i, z)
            k = f_handles.dk_ds(hyp, x, z, i, pdx, pdz) ...
                .* f_handles.ds_dxi(hyp, x, z, i, pdx, pdz);
        else
            % [2.c] k(dx/dx_j, dz/dz_j)
            k = f_handles.d2k_ds2(hyp, x, z, i, pdx, pdz) ...
                .* f_handles.ds_dxi(hyp, x, z, i, pdx, pdz) ...
                .* f_handles.ds_dzj(hyp, x, z, i, pdx, pdz) ...
              + f_handles.dk_ds(hyp, x, z, i, pdx, pdz) ...
                .* f_handles.d2s_dxi_dzj(hyp, x, z, i, pdx, pdz);
        end
    end
    
% derivative w.r.t hyperparameters
else
    if i == 1
        if pdx == 0
            if pdz == 0
                % [1-a] dk(x, z)/dlog(ell)
                k = ell * f_handles.dk_ds(hyp, x, z, i, pdx, pdz) ...
                          .* f_handles.ds_dell(hyp, x, z, i, pdx, pdz);
            else
                % [3.b] dk(x, dz/dz_j)/dlog(ell)
                k = ell * ( f_handles.d2k_ds2(hyp, x, z, i, pdx, pdz) ...
                            .* f_handles.ds_dell(hyp, x, z, i, pdx, pdz) ...
                            .* f_handles.ds_dzj(hyp, x, z, i, pdx, pdz) ...
                          + f_handles.dk_ds(hyp, x, z, i, pdx, pdz) ...
                            .* f_handles.d2s_dell_dzj(hyp, x, z, i, pdx, pdz) );
            end
        else
            if pdz == 0
                % [3.a] dk(dx/dx_i, z)/dlog(ell)
                k = ell * ( f_handles.d2k_ds2(hyp, x, z, i, pdx, pdz) ...
                            .* f_handles.ds_dell(hyp, x, z, i, pdx, pdz) ...
                            .* f_handles.ds_dxi(hyp, x, z, i, pdx, pdz) ...
                          + f_handles.dk_ds(hyp, x, z, i, pdx, pdz) ...
                            .* f_handles.d2s_dell_dxi(hyp, x, z, i, pdx, pdz) );
            else
                % [3.c] dk(dx/dx_j, dz/dz_j)/dlog(ell)
                k = ell * ( f_handles.d3k_ds3(hyp, x, z, i, pdx, pdz) ...
                            .* f_handles.ds_dell(hyp, x, z, i, pdx, pdz) ...
                            .* f_handles.ds_dxi(hyp, x, z, i, pdx, pdz) ...
                            .* f_handles.ds_dzj(hyp, x, z, i, pdx, pdz) ...
                          + f_handles.d2k_ds2(hyp, x, z, i, pdx, pdz) ...
                            .* ( f_handles.d2s_dell_dxi(hyp, x, z, i, pdx, pdz) ...
                                 .* f_handles.ds_dzj(hyp, x, z, i, pdx, pdz) ...
                              + f_handles.ds_dxi(hyp, x, z, i, pdx, pdz) ...
                                .* f_handles.d2s_dell_dzj(hyp, x, z, i, pdx, pdz) ...
                              + f_handles.ds_dell(hyp, x, z, i, pdx, pdz) ...
                                .* f_handles.d2s_dxi_dzj(hyp, x, z, i, pdx, pdz) ) ...
                          + f_handles.dk_ds(hyp, x, z, i, pdx, pdz) ...
                            .* f_handles.d3s_dell_dxi_dzj(hyp, x, z, i, pdx, pdz) );
            end
        end
    else
        % 
        % [1.0] dk(x, z)/dlog(sigma_f)                = 2 * k(x, z)
        % [3.0] dk(x, dz/dz_j)/dlog(sigma_f)          = 2 * k(x, dz/dz_j)
        %       dk(dx/dx_i, z)/dlog(sigma_f)          = 2 * k(dx/dx_i, z)
        %       dk(dx/dx_j, dz/dz_j)/dlog(sigma_f)    = 2  k(dx/dx_j, dz/dz_j)
        k = 2 * covIsoDiff(hyp, x, z, 0, pdx, pdz, f_handles);
    end
end