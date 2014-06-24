function K = covSEisoDerivatives(hyp, x, z, i)

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

ell = exp(hyp(1));                                 % characteristic length scale
sqrt_ell = exp(0.5*hyp(1));
inv_ell2 = exp(-2*hyp(1));
inv_ell4 = exp(-4*hyp(1));

% n, d
[n, d] = size(x);
d = d - 1;

% derivative w.r.t param
deri_wrt_param = 0;
if nargin == 4
    deri_wrt_param = i;
end
        
% covariance matrix
if dg
    K = covSEiso(hyp, x, z);
else
    if xeqz
        K = zeros(n, n);
        
        % from
        for ii = 0:d
            mask1 = x(:, 1) == ii;
            if sum(mask1) > 0
                % to
                for jj = 0:d
                    mask2 = x(:, 1) == jj;
                    if sum(mask2) > 0
                        % pre-calculation
                        K_temp = covSEiso(hyp, x(mask1, 2:end), x(mask2, 2:end));
                        
                        % pair-wise cases
                        if ii == 0                            
                            if jj == 0                                
                                % 1. value-to-value
                                if deri_wrt_param == 1 % derivative w.r.t log ell
                                    K_temp = sq_dist(x(mask1, 2:end)'/ell, x(mask2, 2:end)'/ell).*K_temp;
                                end
                                K(mask1, mask2) = K_temp;
                            else
                                % 2. value-to-deri
                                if deri_wrt_param == 1 % derivative w.r.t log ell
                                    K_temp = (sq_dist(x(mask1, 2:end)'/ell, x(mask2, 2:end)'/ell)-2).*K_temp;
                                end
                                K(mask1, mask2) = inv_ell2*bsxfun(@minus, x(mask1, jj+1), x(mask2, jj+1)').*K_temp;
                            end                            
                        else
                            if jj == 0
                                % 3. deri-to-value
                                if deri_wrt_param == 1 % derivative w.r.t log ell
                                    K_temp = (sq_dist(x(mask1, 2:end)'/ell, x(mask2, 2:end)'/ell)-2).*K_temp;
                                end
                                K(mask1, mask2) = -inv_ell2*bsxfun(@minus, x(mask1, ii+1), x(mask2, ii+1)').*K_temp;
                            else
                                % 4. deri-to-deri
                                if ii == jj
                                    delta = inv_ell2;
                                else
                                    delta = 0;
                                end
                                K_temp2 = (delta - inv_ell4*bsxfun(@minus, x(mask1, ii+1), x(mask2, ii+1)')...
                                                          .*bsxfun(@minus, x(mask1, jj+1), x(mask2, jj+1)')).*K_temp;
                                if deri_wrt_param == 1 % derivative w.r.t log ell
                                    K(mask1, mask2) = (-2*delta + (4*inv_ell4)*bsxfun(@minus, x(mask1, ii+1), x(mask2, ii+1)')...
                                                                             .*bsxfun(@minus, x(mask1, jj+1), x(mask2, jj+1)')).*K_temp ...
                                                      + sq_dist(x(mask1, 2:end)'/ell, x(mask2, 2:end)'/ell).*K_temp2;
                                else
                                    K(mask1, mask2) = K_temp2;
                                end
                            end
                        end
                    end
                end
            end
        end
    else
        ns = size(z, 1);
        K = zeros(n, ns);
        
        % derivatives
        for ii = 0:d
            mask1 = x(:, 1) == ii;
            if sum(mask1) > 0
                K_temp = covSEiso(hyp, x(mask1, 2:end), z);
                        
                % pair-wise cases
                if ii == 0                            
                    % 1. value-to-value
                    K(mask1, :) = K_temp;
                else
                    % 2. deri-to-value
                    K(mask1, :) = -inv_ell2*bsxfun(@minus, x(mask1, ii+1), z(:, ii)').*K_temp;
                end
            end
        end
    end
end

% derivative w.r.t log sf
if deri_wrt_param == 2
    K = 2*K;
end