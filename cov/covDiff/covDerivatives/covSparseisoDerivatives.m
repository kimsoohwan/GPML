function K = covSparseisoDerivatives(hyp, x, z, ii)
%                    2+cos(2*pi*r)         1
% k(x, x') = sf2 * [--------------(1-r) + ---- sin(2*pi*r)]
%                         3               2*pi
%
%     |x-x'|
% r = ------
%      ell

if nargin<2, K = '1'; return; end                  % report number of parameters
% if nargin<2, K = '2'; return; end                  % report number of parameters
if nargin<3, z = []; end                                   % make sure, z exists
xeqz = numel(z)==0; dg = strcmp(z,'diag') && numel(z)>0;        % determine mode
if nargin<4, ii = 0; end                                  % make sure, ii exists

% n, d
[n, d] = size(x);
d = d - 1;

if dg
    %   F*F*
    K = covSparseiso(hyp, x, z);      % FF
else
    if xeqz
        %    |  F  |  D1  |  D2  |  D3  |
        % -------------------------------
        % F  | FF  |  FD1 |  FD2 |  FD3 |
        % D1 |  -  | D1D1 | D1D2 | D1D3 |
        % D2 |  -  |   -  | D2D2 | D2D3 |
        % D3 |  -  |   -  |   -  | D3D3 |
        
        K = zeros(n, n);
        for i = 0:d
            mask_i = x(:, 1) == i;
            for j = i:d
                mask_j = x(:, 1) == j;
                
                % itself
                if i == 0
                    if j == 0
                        K(mask_i, mask_j) = covSparseiso(hyp, x(mask_i, 2:end), x(mask_j, 2:end), ii);      % FF
                    else
                        K(mask_i, mask_j) = covSparseisoFD(hyp, x(mask_i, 2:end), x(mask_j, 2:end), j, ii);	% FD
                    end
                else
                    K(mask_i, mask_j) = covSparseisoDD(hyp, x(mask_i, 2:end), i, x(mask_j, 2:end), j, ii);	% DD
                end
                
                % transpose
                if i ~= j
                    K(mask_j, mask_i) = K(mask_i, mask_j)';
                end
            end
        end
    else
        %    |  F*  |
        % -----------
        % F  | FF*  |
        % D1 | D1F* |
        % D2 | D2F* |
        % D3 | D3F* |
        
        m = size(z, 1);
        K = zeros(n, m);
        for i = 0:d
            mask_i = x(:, 1) == i;
            if i == 0
                K(mask_i, :) = covSparseiso(hyp, x(mask_i, 2:end), z);          % FF
            else
                K(mask_i, :) = covSparseisoFD(hyp, z, x(mask_i, 2:end), i)';    % DF
            end
        end
    end
end