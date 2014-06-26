function K = covSEisoFD(hyp, x, z, j, ii)
% K_SE(f1, \partial f2 / \partial x_j)

% hyperparameters
ell = exp(hyp(1));                              % characteristic length scale
ell2 = exp(2*hyp(1));

% precompute covSEiso
K = covSEiso(hyp, x, z);        % covariance
if ii ~= 0
    dK_theta_i = covSEiso(hyp, x, z, ii);       % derivative of covariance w.r.t hyperparameters 
end

% delta matrix
delta_j = bsxfun(@minus, x(:, j)/ell2, (z(:, j)')/ell2);  % cross delta Kxz

switch ii
    % covariances
    case 0
        K = K .* delta_j;
        
    % derivatives w.r.t log ell
    case 1
        K = delta_j .* (dK_theta_i - 2*K);
        
    % derivatives w.r.t log sf
    case 2
       K = dK_theta_i .* delta_j;
        
    otherwise
      error('Unknown hyperparameter');     
end