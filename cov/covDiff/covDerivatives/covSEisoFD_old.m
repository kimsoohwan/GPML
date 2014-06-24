function K = covSEisoFD(hyp, x, z, j, i)

% Squared Exponential covariance function with isotropic distance measure. The 
% covariance function is parameterized as:
%
% k(x^p,x^q) = sf^2 * exp(-(x^p - x^q)'*inv(P)*(x^p - x^q)/2) 
%
% where the P matrix is ell^2 times the unit matrix and sf^2 is the signal
% variance. The hyperparameters are:
%
% hyp = [ log(ell)
%         log(sf)  ]
%
% For more help on design of covariance functions, try "help covFunctions".
%
% Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch, 2010-09-10.
%
% See also COVFUNCTIONS.M.

if nargin<2, K = '2'; return; end                  % report number of parameters
if nargin<3, z = []; end                                   % make sure, z exists
xeqz = numel(z)==0; dg = strcmp(z,'diag') && numel(z)>0;        % determine mode

ell2 = exp(2*hyp(1));                              % characteristic length scale
    
%
if dg                                                               % vector kxx
  K = zeros(size(x,1),1);
else
    % precompute covSEiso
    K_SE = covSEiso(hyp, x, z);        % covariance
    if nargin == 5
        Partial_K_SE = covSEiso(hyp, x, z, i);     % derivative of covariance w.r.t hyperparameters
    end
    
    % delta matrix
    if xeqz
        K = bsxfun(@minus, x(:, j)/ell2, (x(:, j)')/ell2);  % symmetric matrix Kxx
    else
        K = bsxfun(@minus, x(:, j)/ell2, (z(:, j)')/ell2);  % cross delta Kxz
    end
    
    % covariance matrix       
    if nargin < 5
        K = K .* K_SE;
    else
        K = K .* (Partial_K_SE - 2*K_SE);
    end
end