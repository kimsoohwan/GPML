function K = covisoDerObsBwiseUnstable(f_handles, hyp, x, xd, z, i)

%% function name convention
% cov:      covariance function
% iso:      isotropic
% Diff:     differentiable w.r.t input coordinates or take derivative observations
% Bwise:    block maxtrix-wised version
% Unstable: cannot handle divide-by-zero

%% input/output arguments
% f_handles:    function handles for component functions of derivatives
% hyp:  [1x2]   hyperparameters, hyp = [log(ell), log(sigma_f)]
% x:    [nxd]   first function input vectors
% xd:   [ndxd]  first derivative input vectors
% z:    [nsxd]  second function input vectors, default: [] meaning z = x
% i:            partial deriavtive coordiante w.r.t hyperparameters, default: 0
% K:    [nnxns]  covariance, nn = n + nd*d

% some constants
xeqz = numel(z)==0; dg = strcmp(z,'diag') && numel(z)>0;        % determine mode
[n, d]  = size(x);
ns      = size(z, 1);
nd      = size(xd, 1);

% prediction
if dg
    % Kss
    assert(i == 0);
    K = exp(2*hyp(end))*ones(n+nd, 1); % sigma_f = exp(hyp(end))
    
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

        % number of blocks
        if nd > 0
            num_blocks = 1 + d;
        else
            num_blocks = 1;
        end
        
        % for each row block: x or xd
        for row_block = 1:num_blocks
            % derivative w.r.t x_i
            if row_block == 1
                xx  = x;
                pdx = 0;
                start_row   = 1;
                end_row     = n;
            else
                xx  = xd;
                pdx = row_block-1;
                start_row   = n + nd*(row_block-2) + 1;
                end_row     = n + nd*(row_block-1);
            end
            idx_row = start_row:end_row;

            % for each col: z or zd
            for col_block = 1:num_blocks
                % derivative w.r.t x_i
                if col_block == 1
                    zz  = x;
                    pdz = 0;
                    start_col   = 1;
                    end_col     = n;
                else
                    zz  = xd;
                    pdz = col_block-1;
                    start_col   = n + nd*(col_block-2) + 1;
                    end_col     = n + nd*(col_block-1);
                end
                idx_col = start_col:end_col;

                % calculation
                K(idx_row, idx_col) = covisoDerObsUnstable(f_handles, hyp, xx, zz, i, pdx, pdz);
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
        
        assert(i == 0);
        
        nn  = n  + nd*d;        
        K = zeros(nn, ns); 
        
        % number of blocks
        if nd > 0
            num_blocks = 1 + d;
        else
            num_blocks = 1;
        end
        
        % for each row: x or xd
        for row_block = 1:num_blocks
            % derivative w.r.t x_i
            if row_block == 1
                xx  = x;
                pdx = 0;
                start_row   = 1;
                end_row     = n;
            else
                xx  = xd;
                pdx = row_block-1;
                start_row   = n + nd*(row_block-2) + 1;
                end_row     = n + nd*(row_block-1);
            end
            idx_row = start_row:end_row;

            % for each col: z
            zz = z;
            pdz = 0;
            start_col   = 1;
            end_col     = ns;
            idx_col = start_col:end_col;
            
            % calculation
            K(idx_row, idx_col) = covisoDerObsUnstable(f_handles, hyp, xx, zz, i, pdx, pdz);
        end
    end
end
