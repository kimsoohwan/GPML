function K = covSEisoFD(hyp, x, z, j, ii)
% K_SE(f1, \partial f2 / \partial x_j)

if nargin<2, K = '2'; return; end                  % report number of parameters
if nargin<5, ii = 0; end                                      % derivative index

% hyperparameters
ell2 = exp(2*hyp(1));                              % characteristic length scale

% precompute covSEiso
K_SE = covSEiso(hyp, x, z);        % covariance
if nargin == 5
    if ii == 0,    Partial_K_SE = covSEiso(hyp, x, z);              % covariance
    else           Partial_K_SE = covSEiso(hyp, x, z, ii); end      % derivative of covariance w.r.t hyperparameters 
end

% delta matrix
K = bsxfun(@minus, x(:, j)/ell2, (z(:, j)')/ell2);  % cross delta Kxz

switch ii
    % covariances
    case 0
        K = K .* K_SE;
        
    % derivatives w.r.t log ell
    case 1
        K = K .* (Partial_K_SE - 2*K_SE);
        
    % derivatives w.r.t log sf
    case 2
       %K =     K .* Partial_K_SE; Bug!!!
        K = 2 * K .* Partial_K_SE;
        
    otherwise
      error('Unknown hyperparameter');     
end