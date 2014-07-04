function K = covSparseiso(hyp, x, z, ii)

%                    2+cos(2*pi*r)         1
% k(x, x') = sf2 * [--------------(1-r) + ---- sin(2*pi*r)]
%                         3               2*pi
%
%     |x-x'|
% r = ------
%      ell
%
% [X,Y] = meshgrid(-2:.05:2, -2:.05:2);                                
% R = sqrt(X.^2 + Y.^2);      
% IV = R >= 1;                                  
% Z = (2+cos(2*pi*R)).*(1-R)/3 + sin(2*pi*R)/(2*pi);
% Z(IV) = 0;
% surf(X,Y,Z)

if nargin<2, K = '2'; return; end                  % report number of parameters
if nargin<3, z = []; end                                   % make sure, z exists
xeqz = numel(z)==0; dg = strcmp(z,'diag') && numel(z)>0;        % determine mode
if nargin<4, ii = 0; end                                  % make sure, ii exists

ell = exp(hyp(1));
if length(hyp) == 1
    sf2 = 1;
else
    sf2 = exp(2*hyp(2));
end

% precompute distances
% R = sqrt((x-z)'(x-z)) / ell = r/ell
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

switch ii
    case 0  % covariances
        K = sf2*((2+cos(2*pi*R)).*(1-R)/3 + sin(2*pi*R)/(2*pi));
    case 1  % derivatives w.r.t log ell
        pK_pR = (2*sf2/3)*(cos(2*pi*R) - pi*sin(2*pi*R).*(1-R)-1);
        pR_pL = -R/ell;
        K = ell*pK_pR.*pR_pL;
        %K = -(2*sf2/3)*(cos(2*pi*R) - pi*sin(2*pi*R).*(1-R)-1).*R;
    case 2  % derivatives w.r.t log sf 
        K = 2*sf2*((2+cos(2*pi*R)).*(1-R)/3 + sin(2*pi*R)/(2*pi));
    otherwise
        error('Unknown hyperparameter')
end

% remove invalid coefficients
K(invalid_mask) = 0;