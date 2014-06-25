function K = covSEisoDD(hyp, x, i, z, j, ii)
% K_SE(\partial f1 / \partial x_i, \partial f2 / \partial x_j)

% hyperparameters
ell2 = exp(2*hyp(1));                              % characteristic length scale

% precompute covSEiso
K = covSEiso(hyp, x, z);        % covariance
if ii ~= 0
    dK_theta_i = covSEiso(hyp, x, z, ii);     % derivative of covariance w.r.t hyperparameters
end

% delta(i, j)
if i == j
    delta = 1;
else
    delta = 0;
end

% delta matrix
delta_i = bsxfun(@minus, x(:, i)/ell2, (z(:, i)')/ell2);  % cross delta Kxz
delta_j = bsxfun(@minus, x(:, j)/ell2, (z(:, j)')/ell2);  % cross delta Kxz

% covariance matrix
factor = delta/ell2 - delta_i.*delta_j; % factor

switch ii
    % covariances
    case 0
        K = K .* factor;
      
    % derivatives w.r.t log ell
    case 1
        K = dK_theta_i .* factor + K .* (-2*delta/ell2 + 4*delta_i.*delta_j);
      
    % derivatives w.r.t log sf
    case 2
        K = dK_theta_i .* factor;
        
    otherwise
        error('Unknown hyperparameter')
end