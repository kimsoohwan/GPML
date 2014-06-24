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
% R = sqrt((x-x')'(x-x')) / ell
if dg                                                               % vector kxx
  R = zeros(size(x,1),1);
else
  if xeqz                                                 % symmetric matrix Kxx
    R = sqrt( sq_dist(x'/ell) );
  else                                                   % cross covariances Kxz
    R = sqrt( sq_dist(x'/ell, z'/ell) );
  end
end

% invalid mask
invalid_mask = R >= 1; % r >= 1

% avoiding division by 0
mask_R0 = R < eps; 

% \partial k / \partial r
pK_pR	= (2*sf2/3)*(cos(2*pi*R) - pi*sin(2*pi*R).*(1-R)-1);
p2K_pR2	= (-2*pi*sf2/3)*(sin(2*pi*R) + 2*pi*cos(2*pi*R).*(1-R));
p3K_pR3	= (8*pi^3*sf2/3)*sin(2*pi*R).*(1-R);

% delta matrix:
pR_pi =  bsxfun(@minus, x(:, i)/ell2, (z(:, i)')/ell2)./R;  pR_pi(mask_R0) = 0;
pR_pj = -bsxfun(@minus, x(:, j)/ell2, (z(:, j)')/ell2)./R;  pR_pj(mask_R0) = 0;

% delta(i, j)
if i == j
    delta = 1;
else
    delta = 0;
end
p2R_pipj    = - delta./(ell2*R) - pR_pi.*pR_pj./R; p2R_pipj(mask_R0) = 0;

pR_pL       = -R/ell;
p2R_pipL    = -pR_pi/ell;
p2R_pjpL    = -pR_pj/ell;
p3R_pipjpL  = -p2R_pipj/ell;

switch ii
    case 0  % covariances
        K = p2K_pR2.*pR_pi.*pR_pj + pK_pR.*p2R_pipj;
    case 1  % derivatives w.r.t log ell
        K = ell*(p3K_pR3.*pR_pi.*pR_pj.*pR_pL + p2K_pR2.*(p2R_pipj.*pR_pL + pR_pi.*p2R_pjpL + pR_pj.*p2R_pipL) + pK_pR.*p3R_pipjpL);
        %K = -p2R_pipj.*(pK_pR + R.*p2K_pR2) - pR_pi.*pR_pj.*(2*p2K_pR2 + R.*p3K_pR3);
    case 2  % derivatives w.r.t log sf 
        K = 2*(p2K_pR2.*pR_pi.*pR_pj + pK_pR.*p2R_pipj);
    otherwise
        error('Unknown hyperparameter')
end

% remove invalid coefficients
K(invalid_mask) = 0;