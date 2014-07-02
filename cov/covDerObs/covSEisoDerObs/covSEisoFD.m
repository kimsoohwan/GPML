function K = covSEisoFD(hyp, x, z, j, ii)
% K_SE(f1, \partial f2 / \partial x_j)

% hyperparameters
ell = exp(hyp(1));                              % characteristic length scale
ell2 = exp(2*hyp(1));

% covSEiso
K = covSEiso(hyp, x, z);        % covariance

% delta matrix
delta_j = bsxfun(@minus, x(:, j)/ell2, (z(:, j)')/ell2);  % cross delta Kxz

% K
K = K .* delta_j;

switch ii
    % covariances
    case 0

    % derivatives w.r.t log ell
    case 1
        S = (-1/2)*sq_dist(x'/ell, z'/ell);
        K = K .* (-2*S - 2);
        
    % derivatives w.r.t log sf
    case 2
       K = 2*K;
        
    otherwise
      error('Unknown hyperparameter');     
end