function K = covDiff(hyp, x, xd, z, ii, hCovFF, hCovFD, hCovDD)

if nargin<4, z = []; end                                   % make sure, z exists
xeqz = numel(z)==0; dg = strcmp(z,'diag') && numel(z)>0;        % determine mode
if nargin<5, ii = 0; end                                      % derivative index

% number of data and dimensions
[n, d] = size(x);
n1 = n - n2;
x1 = x(1:n1, :);
x2 = x(n1+1:end, :);

if dg
    K = feval(hCovFF{:}, hyp, x, z);
else
    if xeqz
        %            |F1(n1)|D1(n1)|D2(n1)|D3(n1)| F2(n2)
        % K = ------------------------------------------
        %     F1(n1) | F1F1,  F1D1,  F1D2,  F1D3,| F1F2
        %     D1(n1) |    -,  D1D1,  D1D2,  D1D3,| D1F2
        %     D2(n1) |    -,     -,  D2D2,  D2D3,| D2F2
        %     D3(n1) |    -,     -,     -,  D3D3,| D3F2
        %     ------------------------------------------
        %     F2(n2) |    -,     -,     -,    -, | F2F2
        
        K = zeros(n1*(d+1)+n2, n1*(d+1)+n2);
        for row = 0:d
            idx_row = n1*row+1:n1*(row+1);
            for col = row:d
                idx_col = n1*col+1:n1*(col+1);
                
                % itself
                if row == 0
                    if col == 0
                        % F1F1
                        if ii == 0,	K(idx_row, idx_col) = feval(hCovFF{:}, hyp, x1);
                        else        K(idx_row, idx_col) = feval(hCovFF{:}, hyp, x1, [], ii); end
                    else
                        % F1D*
                        K(idx_row, idx_col) = feval(hCovFD{:}, hyp, x1, x1, col, ii);
                    end
                else
                        % D*D*
                        K(idx_row, idx_col) = feval(hCovDD{:}, hyp, x1, row, x1, col, ii);
                end

                % transpose
                if row ~= col
                    K(idx_col, idx_row) = K(idx_row, idx_col)';
                end
                
            end
            % last column
            idx_col = n1*(d+1)+1:n1*(d+1)+n2;
            
            % F1F2
            if row == 0
                if ii == 0,	K(idx_row, idx_col) = feval(hCovFF{:}, hyp, x1, x2);
                else        K(idx_row, idx_col) = feval(hCovFF{:}, hyp, x1, x2, ii); end
                
            % D*F2
            else
                K(idx_row, idx_col) = feval(hCovFD{:}, hyp, x2, x1, row, ii)';
            end
            
            % transpose
            K(idx_col, idx_row) = K(idx_row, idx_col)';
        end
        % last row and last column
        idx_row = n1*(d+1)+1:n1*(d+1)+n2;
        idx_col = n1*(d+1)+1:n1*(d+1)+n2;
        
        % F2F2
        if ii == 0,	K(idx_row, idx_col) = feval(hCovFF{:}, hyp, x2, []);
        else        K(idx_row, idx_col) = feval(hCovFF{:}, hyp, x2, [], ii); end        
    else
        %            |F(m)|
        % K = -------------
        %     F1(n1) | F1F
        %     D1(n1) | D1F
        %     D2(n1) | D2F
        %     D3(n1) | D3F
        %     -------------
        %     F2(n2) | F2F
        
        [m, d] = size(z);
        K = zeros(n1*(d+1)+n2, m);
        for row = 0:d
            idx_row = n1*row+1:n1*(row+1);
            if row == 0
                % F1F
                if ii == 0,	K(idx_row, :) = feval(hCovFF{:}, hyp, x1, z);
                else        K(idx_row, :) = feval(hCovFF{:}, hyp, x1, z, ii); end
            else
                % D*F
                K(idx_row, :) = feval(hCovFD{:}, hyp, z, x1, row, ii)';
            end
        end
        % last row
        idx_row = n1*(d+1)+1:n1*(d+1)+n2;
        
        % F2F
        if ii == 0,	K(idx_row, :) = feval(hCovFF{:}, hyp, x2, z);
        else        K(idx_row, :) = feval(hCovFF{:}, hyp, x2, z, ii); end
    end
end