function K = covMaterniso3FD(hyp, x, z, j, ii)
% K_Materniso3(f1, \partial f2 / \partial x_j)

if nargin<2, K = '2'; return; end                  % report number of parameters
if nargin<5, ii = 0; end                                      % derivative index

% hyperparameters
ell  = exp(hyp(1));                              % characteristic length scale
ell2 = exp(2*hyp(1));
sf2  = exp(2*hyp(2));

% precompute S = -sqrt(3)*r/ell
S = -sqrt(3*sq_dist(x'/ell, z'/ell));

% delta matrix: (x_j - x'_j)/ell^2
delta_j = bsxfun(@minus, x(:, j)/ell2, (z(:, j)')/ell2);  % cross delta Kxz

% k = sf2 * (1 + sqrt(3)*r/ell) * exp(-sqrt(3)*r/ell)
%   = sf2 * (1 - s) * exp(s)
%
% s^2 = (3/ell^2) * sum_{i=1}^d (xi - zi)^2
% ds/dzj = -(3/ell^2)*(xj - zj)/s
% dk/dzj = sf2 * exp(s) * (-ds/dzj + (1-s)*ds/dzj)
%        = sf2 * exp(s) * (-s*ds/dzj)
%        = 3*sf2 * exp(s) * (xj - zj)/ell^2
K = (3*sf2) * exp(S) .* delta_j;

switch ii
    % covariances
    case 0

    % derivatives w.r.t log ell
    case 1
        % ds/dell = (-1/ell)*s
        % d2k/dlog(ell) dzj = ell * 3*sf2 * exp(s) * (xj-zj)/ell^2 * (-2/ell -s/ell)
        %                   = 3*sf2 * exp(s) * (xj-zj)/ell^2 * (-s -2)
        K = K.*(-S - 2);
        
    % derivatives w.r.t log sf
    case 2
        K = 2*K;
        
    otherwise
      error('Unknown hyperparameter');     
end