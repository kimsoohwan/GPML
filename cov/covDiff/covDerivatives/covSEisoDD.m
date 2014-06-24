function K = covSEisoDD(hyp, x, i, z, j, ii)
% K_SE(\partial f1 / \partial x_i, \partial f2 / \partial x_j)

if nargin<2, K = '2'; return; end                  % report number of parameters
if nargin<6, ii = 0; end                                      % derivative index

% hyperparameters
ell2 = exp(2*hyp(1));                              % characteristic length scale

% precompute covSEiso
K_SE = covSEiso(hyp, x, z);        % covariance
if ii ~= 0
    Partial_K_SE = covSEiso(hyp, x, z, ii);     % derivative of covariance w.r.t hyperparameters
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
K = delta/ell2 - delta_i.*delta_j; % factor

switch ii
    % covariances
    case 0
        K = K.*K_SE;
      
    % derivatives w.r.t log ell
    case 1
        K = (-2*delta/ell2 + 4*delta_i.*delta_j).*K_SE + K.*Partial_K_SE;
      
    % derivatives w.r.t log sf
    case 2
        K = K.*Partial_K_SE;
        
    otherwise
        error('Unknown hyperparameter')
end