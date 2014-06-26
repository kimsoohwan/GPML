function K = covisoDiffCwiseUnstable(f_handles, hyp, x, xd, z, i)

%% function name convention
% cov:  covariance function
% iso:  isotropic
% Diff: differentiable w.r.t input coordinates or take derivative observations
% Cwise: non-vectorized or coefficient-wised version
% Unstable: cannot handle divide-by-zero

%% input/output arguments
% hyp:  [1x2]   hyperparameters, hyp = [log(ell), log(sigma_f)]
% x:    [nxd]   first function input vectors
% z:    [nsxd]  second function input vectors, default: [] meaning z = x
% i:            partial deriavtive coordiante w.r.t hyperparameters, default: 0
% xd:   [ndxd]  first derivative input vectors
% K:    [nnxns]  covariance, nn = n + nd*d

% some constants
xeqz = numel(z)==0; dg = strcmp(z,'diag') && numel(z)>0;        % determine mode
[n, d]  = size(x);
ns      = size(z, 1);
nd      = size(xd, 1);

% prediction
if dg
    % Kss
    K = ones(n+nd, 1);
    
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
                K(row, col) = covisoDiffUnstable(f_handles, hyp, xx, zz, i, pdx, pdz);
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
                K(row, col) = covisoDiffUnstable(f_handles, hyp, xx, zz, i, pdx, pdz);
            end
        end
    end
end
