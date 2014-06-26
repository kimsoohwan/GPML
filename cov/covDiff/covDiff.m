function K = covDiff(f_handles, hyp, x, xd, z, i)

%% function name convention
% cov:  covariance function
% Diff: differentiable w.r.t input coordinates or take derivative observations
% Bwise: block maxtrix-wised version

%% input/output arguments
% f_handles:    function handles for component functions of derivatives
% hyp:          hyperparameters
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
        % f(x)        : n  |  FF  |  FD1,  FD2,  FD3
        % df(xd)/dx_1 : nd |   -  | D1D1, D1D2, D1D3  
        % ...              |
        % df(xd)/dx_d : nd |   -  |   - ,  ..., D3D3
        
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
                start_row   = 1;
                end_row     = n;
            else
                start_row   = n + nd*(row_block-2) + 1;
                end_row     = n + nd*(row_block-1);
            end
            idx_row = start_row:end_row;

            % for each col: z or zd
            for col_block = row_block:num_blocks
                % derivative w.r.t x_i
                if col_block == 1
                    start_col   = 1;
                    end_col     = n;
                else
                    start_col   = n + nd*(col_block-2) + 1;
                    end_col     = n + nd*(col_block-1);
                end
                idx_col = start_col:end_col;
                
                % calculation
                
                % itself
                if row_block == 1
                    if col_block == 1
                        % FF
                        if i == 0,	K(idx_row, idx_col) = feval(f_handles.FF{:}, hyp, x);
                        else        K(idx_row, idx_col) = feval(f_handles.FF{:}, hyp, x, [], i); end
                    else
                        % FD*
                        K(idx_row, idx_col) = feval(f_handles.FD{:}, hyp, x, xd, col_block-1, i);                        
                    end
                else
                        % D*D*
                        K(idx_row, idx_col) = feval(f_handles.DD{:}, hyp, xd, row_block-1, xd, col_block-1, i);
                end
                
                % transpose
                if row_block ~= col_block
                    K(idx_col, idx_row) = K(idx_row, idx_col)';
                end
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
                start_row   = 1;
                end_row     = n;
            else
                start_row   = n + nd*(row_block-2) + 1;
                end_row     = n + nd*(row_block-1);
            end
            idx_row = start_row:end_row;

            % for each col: z
            start_col   = 1;
            end_col     = ns;
            idx_col = start_col:end_col;
            
            % calculation
            if row_block == 1
                % FF
                K(idx_row, idx_col) = feval(f_handles.FF{:}, hyp, x, z);
            else
                % D*F
                K(idx_row, idx_col) = feval(f_handles.FD{:}, hyp, z, xd, row_block-1, i)';
            end
        end
    end
end
