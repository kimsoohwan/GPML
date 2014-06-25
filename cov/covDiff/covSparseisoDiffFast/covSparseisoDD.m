function K = covSparseisoDD(hyp, x, i, z, j, ii)

%                    2+cos(2*pi*r)         1
% k(x, x') = sf2 * [--------------(1-r) + ---- sin(2*pi*r)]
%                         3               2*pi
%
%     |x-x'|
% r = ------
%      ell

if nargin<2, K = '2'; return; end                  % report number of parameters
if nargin<4, z = []; end                                   % make sure, z exists
xeqz = numel(z)==0; dg = strcmp(z,'diag') && numel(z)>0;        % determine mode
if nargin<6, ii = 0; end                                  % make sure, ii exists

ell = exp(hyp(1));
ell2 = exp(2*hyp(1));
if length(hyp) == 1
    sf2 = 1;
else
    sf2 = exp(2*hyp(2));
end

% precompute distances
% S = sqrt((x-x')'(x-x')) / ell
if dg                                                               % vector kxx
  S = zeros(size(x,1),1);
else
  if xeqz                                                 % symmetric matrix Kxx
    S = sqrt( sq_dist(x'/ell) );
  else                                                   % cross covariances Kxz
    S = sqrt( sq_dist(x'/ell, z'/ell) );
  end
end

% invalid mask
invalid_mask = S >= 1; % r >= 1

% avoiding division by 0
mask_S0 = S < eps; 

% \partial k / \partial r
pK_pS	= (2*sf2/3)*(cos(2*pi*S) - pi*sin(2*pi*S).*(1-S)-1);
p2K_pS2	= (-2*pi*sf2/3)*(sin(2*pi*S) + 2*pi*cos(2*pi*S).*(1-S));
p3K_pS3	= (8*pi^3*sf2/3)*sin(2*pi*S).*(1-S);

% delta matrix:
delta_i = bsxfun(@minus, x(:, i)/ell2, (z(:, i)')/ell2);
delta_j = bsxfun(@minus, x(:, j)/ell2, (z(:, j)')/ell2);

% delta(i, j)
if i == j
    delta = 1/ell^2;
else
    delta = 0;
end

% K
K = - p2K_pS2 .* delta_i .* delta_j ./ (S.^2) ...
    + pK_pS .* (-delta./S + delta_i .* delta_j ./ (S.^3));
K(mask_S0) = (4/3)*pi^2*sf2*delta;
switch ii
    case 0  % covariances
        
    case 1  % derivatives w.r.t log ell
        K = p3K_pS3 .* delta_i .* delta_j ./ S ...
          + p2K_pS2 .* (delta      + delta_i .* delta_j ./ (S.^2)) ...
          + pK_pS   .* (delta ./ S - delta_i .* delta_j ./ (S.^3));
        K(mask_S0) = -(8/3)*(pi^2)*sf2 * delta;
        
    case 2  % derivatives w.r.t log sf 
        K = 2*K;
    otherwise
        error('Unknown hyperparameter')
end

% remove invalid coefficients
K(invalid_mask) = 0;