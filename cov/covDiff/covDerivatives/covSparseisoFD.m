function K = covSparseisoFD(hyp, x, z, j, ii)

%                    2+cos(2*pi*r)         1
% k(x, x') = sf2 * [--------------(1-r) + ---- sin(2*pi*r)]
%                         3               2*pi
%
%     |x-x'|
% r = ------
%      ell

% [X,Y] = meshgrid(-2:.05:2, -2:.05:2);                                
% R = sqrt(X.^2 + Y.^2);                                        
% Z = (2/3)*(cos(2*pi*R) - pi*sin(2*pi*R).*(1-R)-1).*(-X)./R;
% Z(R<eps) = 0;
% surf(X,Y,Z)
% xlabel('x'); ylabel('y');
%
% [X,Y] = meshgrid(-2:.05:2, -2:.05:2);                                
% R = sqrt(X.^2 + Y.^2); 
% pK_pr = (2/3)*(cos(2*pi*R) - pi*sin(2*pi*R).*(1-R)-1);
% ppK_prpr = (-2*pi/3)*(sin(2*pi*R) + 2*pi*cos(2*pi*R).*(1-R));                                 
% pR_pj = -X./R;
% Z = pR_pj.*(pK_pr + R.*ppK_prpr);
% Z(R<eps) = 0;
% surf(X,Y,Z)
% xlabel('x'); ylabel('y');

if nargin<2, K = '2'; return; end                  % report number of parameters
if nargin<5, ii = 0; end                                  % make sure, ii exists

ell = exp(hyp(1));
ell2 = exp(2*hyp(1));
if length(hyp) == 1
    sf2 = 1;
else
    sf2 = exp(2*hyp(2));
end

% precompute distances
% R = sqrt((x-x')'(x-x')) / ell
R = sqrt( sq_dist(x'/ell, z'/ell) );    % cross covariances Kxz

% invalid mask
invalid_mask = R >= 1; % r >= 1

% avoiding division by 0
mask_R0 = R < eps; 

% \partial k / \partial r
pK_pR   = (2*sf2/3)*(cos(2*pi*R) - pi*sin(2*pi*R).*(1-R)-1);
p2K_pR2 = (-2*pi*sf2/3)*(sin(2*pi*R) + 2*pi*cos(2*pi*R).*(1-R));

% delta matrix: -(x_j - x'_j)/(l^2r)
pR_pj = -bsxfun(@minus, x(:, j)/ell2, (z(:, j)')/ell2)./R;  pR_pj(mask_R0) = 0;
pR_pL       = -R/ell;
p2R_pjpL    = -pR_pj/ell;

%K = pK_pR.*pR_pj;
switch ii
    case 0  % covariances
        K = pK_pR.*pR_pj;
    case 1  % derivatives w.r.t log ell
        K = ell*(p2K_pR2.*pR_pj.*pR_pL + pK_pR.*p2R_pjpL);
        %K = -R.*p2K_pR2.*pR_pj - K;
    case 2  % derivatives w.r.t log sf 
        %K = 2*pr_pj.*pK_pr;
        K = 2*pK_pR.*pR_pj;
    otherwise
        error('Unknown hyperparameter')
end

% remove invalid coefficients
K(invalid_mask) = 0;