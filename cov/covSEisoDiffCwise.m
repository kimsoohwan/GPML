function K = covSEisoDiffCwise(hyp, x, z, i, xd)

%% function name convention
% cov:  covariance function
% SE:   squared exponential
% iso:  isotropic
% Diff: differentiable w.r.t input coordinates or take derivative observations
% Cwise: non-vectorized or coefficient-wised version

%% input/output arguments
% hyp:  [1x2]   hyperparameters, hyp = [log(ell), log(sigma_f)]
% x:    [nxd]   first function/derivative input vector
% z:    [nsxd]  second function/derivative input vector, default: [] meaning z = x
% i:            partial deriavtive coordiante w.r.t hyperparameters, default: 0
% pdx:          partial derivative coordinate w.r.t the first input
%               if dxi = 0 (default), x = function input vector, else x = derivative input
% pdz:          partial derivative coordinate w.r.t the second input
%               if pdz = 0 (default), z = function input vector, else z = derivative input
% K:    [nxns]  covariance

% default parameters
if nargin < 2, K = '2'; return; end                  % report number of parameters
if nargin < 3, z = [];  end                                   % make sure, z exists
if nargin < 4, i = 0;   end
if nargin < 5, xd = []; end
xeqz = numel(z)==0; dg = strcmp(z,'diag') && numel(z)>0;        % determine mode
[n, d]  = size(x);
ns      = size(z, 1);
nd      = size(xd, 1);

% prediction
if dg
    % Kss
    K = zeros(n+nd, 1);
    
% learning
else
    % K
    if xeqz  
        % K
        %                  | f(x) | df(xd)/dx_1, ..., df(xd)/dx_d
        %                  |  n   |    nd                  nd
        % -------------------------------------------------------
        % f(x)        : n  |
        % df(xd)/dx_1 : nd |
        % ...              |
        % df(xd)/dx_d : nd |
        
        nn = n + nd*d;        
        K = zeros(nn, nn);

        % for each row: x or xd
        for row = 1:nn
            % derivative w.r.t x_i
            if row <= n
                xx  = x(row, :);
                pdx = 0;
            else
                xx_index = mod(row-n, nd);
                if xx_index == 0
                    xx_index = nd;
                end
                xx  = xd(xx_index, :);
                pdx = ceil((row-n)/nd);
            end

            % for each col: z or zd
            for col = 1:nn
                % derivative w.r.t x_i
                if col <= n
                    zz  = x(col, :);
                    pdz = 0;
                else
                    zz_index = mod(col-n, nd);
                    if zz_index == 0
                        zz_index = nd;
                    end
                    zz  = xd(zz_index, :);
                    pdz = ceil((col-n)/nd);
                end

                % calculation
                K(row, col) = covSEisoDiffBase(hyp, xx, zz, i, pdx, pdz);
            end
        end
        
    % Ks
    else        
        % K
        %             | f(z)
        % -------------------
        % f(x)        |
        % df(xd)/dx_1 |
        % ...         |
        % df(xd)/dx_d |
        
        nn  = n  + nd*d;        
        K = zeros(nn, ns); 
        
        % for each row: x or xd
        for row = 1:nn
            % derivative w.r.t x_i
            if row <= n
                xx  = x(row, :);
                pdx = 0;
            else
                xx_index = mod(row-n, nd);
                if xx_index == 0
                    xx_index = nd;
                end
                xx  = xd(xx_index, :);
                pdx = ceil((row-n)/nd);
            end

            % for each col: z
            pdz = 0;
            for col = 1:ns
                zz  = z(col, :);

                % calculation
                K(row, col) = covSEisoDiffBase(hyp, xx, zz, i, pdx, pdz);
            end
        end
    end
end
