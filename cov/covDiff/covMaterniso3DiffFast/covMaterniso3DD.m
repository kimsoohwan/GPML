function K = covMaterniso3DD(hyp, x, i, z, j, ii)
% K_Materniso3(\partial f1 / \partial x_i, \partial f2 / \partial x_j)

if nargin<2, K = '2'; return; end                  % report number of parameters
if nargin<6, ii = 0; end                                      % derivative index

% hyperparameters
ell  = exp(hyp(1));                              % characteristic length scale
ell2 = exp(2*hyp(1));
sf2  = exp(2*hyp(2));

% precompute R = sqrt(3)*r/ell
R = sqrt(3*sq_dist(x'/ell, z'/ell));
mask_R0 = R < eps; % avoiding division by 0

% delta(i, j)
if i == j
    delta = 1;
else
    delta = 0;
end

% delta matrix
delta_i = bsxfun(@minus, x(:, i)/ell2, (z(:, i)')/ell2);  % cross delta Kxz
delta_j = bsxfun(@minus, x(:, j)/ell2, (z(:, j)')/ell2);  % cross delta Kxz

% covariances
K = (3*sf2) * exp(-R) .* (delta/ell2 - (3./R).*delta_i.*delta_j);
switch ii
    % covariances
    case 0
        %K = 3*sf2(delta/ell2 - 3*delta_i.*delta_j./R).*exp(-R);
        K(mask_R0) = 3*sf2*delta/ell2;
      
    % derivatives w.r.t log ell
    case 1
        K = (3*sf2) * exp(-R) .* (-2*delta/ell2 + (9./R).*delta_i.*delta_j) + K.*R;
        K(mask_R0) = -6*sf2*delta/ell2;
      
    % derivatives w.r.t log sf
    case 2
        K = 2*K;
        K(mask_R0) = 6*sf2*delta/ell2;
        
    otherwise
        error('Unknown hyperparameter')
end