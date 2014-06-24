function K = covMaterniso3FD(hyp, x, z, j, ii)
% K_Materniso3(f1, \partial f2 / \partial x_j)

if nargin<2, K = '2'; return; end                  % report number of parameters
if nargin<5, ii = 0; end                                      % derivative index

% hyperparameters
ell  = exp(hyp(1));                              % characteristic length scale
ell2 = exp(2*hyp(1));
sf2  = exp(2*hyp(2));

% precompute R = sqrt(3)*r/ell
R = sqrt(3*sq_dist(x'/ell, z'/ell));

% delta matrix: (x_j - x'_j)/l^2
K = bsxfun(@minus, x(:, j)/ell2, (z(:, j)')/ell2);  % cross delta Kxz

K = 3*sf2*K.*exp(-R);
switch ii
    % covariances
    case 0
        %K = 3*sf2*K.*exp(-R);
        
    % derivatives w.r.t log ell
    case 1
        %K = 3*sf2*K.*(R - 2).*exp(-R);
        K = K.*(R - 2);
        
    % derivatives w.r.t log sf
    case 2
        %K = 6*sf2*K.*exp(-R);
        K = 2*K;
        
    otherwise
      error('Unknown hyperparameter');     
end