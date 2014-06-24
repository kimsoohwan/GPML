function K = covMaternisoDerivatives(hyp, x, z, i)
% d = 3

% Matern covariance function with nu = d/2 and isotropic distance measure. For
% d=1 the function is also known as the exponential covariance function or the 
% Ornstein-Uhlenbeck covariance in 1d. The covariance function is:
%
%   k(x^p,x^q) = s2f * f( sqrt(d)*r ) * exp(-sqrt(d)*r)
%
% with f(t)=1 for d=1, f(t)=1+t for d=3 and f(t)=1+t+t²/3 for d=5.
% Here r is the distance sqrt((x^p-x^q)'*inv(P)*(x^p-x^q)), P is ell times
% the unit matrix and sf2 is the signal variance. The hyperparameters are:
%
% hyp = [ log(ell)
%         log(sqrt(sf2)) ]
%
% Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch, 2010-09-10.
%
% See also COVFUNCTIONS.M.


if nargin<2, K = '2'; return; end                  % report number of parameters
if nargin<3, z = []; end                                   % make sure, z exists
xeqz = numel(z)==0; dg = strcmp(z,'diag') && numel(z)>0;        % determine mode

ell = exp(hyp(1));                                 % characteristic length scale
sf2 = exp(2*hyp(2));                                              % signal noise
inv_ell2 = exp(-2*hyp(1));
sqrt3 = sqrt(3);
sqrt3_inv_ell = sqrt3/ell;
three_sf2_inv_ell2 = 3*sf2*inv_ell2;

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
    K = covMaterniso(3, hyp, x, z);
else
    if xeqz
        K = zeros(n, n);
        
        % derivatives
        for ii = 0:d
            mask1 = x(:, 1) == ii;
            if sum(mask1) > 0
                % derivatives
                for jj = 0:d
                    mask2 = x(:, 1) == jj;
                    if sum(mask2) > 0
                        sqrt3_r_inv_l = sqrt(sq_dist(sqrt3_inv_ell*x(mask1, 2:end)', sqrt3_inv_ell*x(mask2, 2:end)')); % sqrt{3}r/l
                        
                        % pair-wise cases
                        if ii == 0                            
                            if jj == 0                                
                                % 1. value-to-value
                                if deri_wrt_param == 1 % derivative w.r.t log ell
                                    K(mask1, mask2) = sf2*(sqrt3_r_inv_l.^2).*exp(-sqrt3_r_inv_l);
                                else
                                    K(mask1, mask2) = sf2*(1 + sqrt3_r_inv_l).*exp(-sqrt3_r_inv_l);
                                end
                            else
                                % 2. value-to-deri
                                if deri_wrt_param == 1 % derivative w.r.t log ell
                                    K(mask1, mask2) = three_sf2_inv_ell2*bsxfun(@minus, x(mask1, jj+1), x(mask2, jj+1)').*(sqrt3_r_inv_l - 2).*exp(-sqrt3_r_inv_l);
                                else
                                    K(mask1, mask2) = three_sf2_inv_ell2*bsxfun(@minus, x(mask1, jj+1), x(mask2, jj+1)').*exp(-sqrt3_r_inv_l);
                                end
                            end                            
                        else
                            if jj == 0
                                % 3. deri-to-value
                                if deri_wrt_param == 1 % derivative w.r.t log ell
                                    K(mask1, mask2) = (-three_sf2_inv_ell2)*bsxfun(@minus, x(mask1, ii+1), x(mask2, ii+1)').*(sqrt3_r_inv_l - 2).*exp(-sqrt3_r_inv_l);
                                else
                                    K(mask1, mask2) = (-three_sf2_inv_ell2)*bsxfun(@minus, x(mask1, ii+1), x(mask2, ii+1)').*exp(-sqrt3_r_inv_l);
                                end
                            else
                                % 4. deri-to-deri
                                if ii == jj
                                    delta = 1;
                                else
                                    delta = 0;
                                end
                                
                                % what if r == 0? -> divide by 0
                                sqrt3_r_inv_l(sqrt3_r_inv_l < eps) = eps;
                                sqrt3_inv_lr = (3*inv_ell2)./sqrt3_r_inv_l; % sqrt(3)/(lr) = (3/l^2) * (l/(sqrt{3}r))
                                
                                if deri_wrt_param == 1 % derivative w.r.t log ell
                                    K(mask1, mask2) = three_sf2_inv_ell2*(((sqrt3_r_inv_l - 2)*delta - sqrt3_inv_lr.*bsxfun(@minus, x(mask1, ii+1), x(mask2, ii+1)') ...
                                                                                                                   .*bsxfun(@minus, x(mask1, jj+1), x(mask2, jj+1)').*(sqrt3_r_inv_l - 3)) ...
                                                                         .*exp(-sqrt3_r_inv_l));
                                else
                                    K(mask1, mask2) = three_sf2_inv_ell2*((delta - sqrt3_inv_lr.*bsxfun(@minus, x(mask1, ii+1), x(mask2, ii+1)') ...
                                                                                               .*bsxfun(@minus, x(mask1, jj+1), x(mask2, jj+1)')) ...
                                                                          .*exp(-sqrt3_r_inv_l));
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
                sqrt3_r_inv_l = sqrt(sq_dist(sqrt3*x(mask1, 2:end)'/ell, sqrt3*z'/ell)); % sqrt{3}r/l
                        
                % pair-wise cases
                if ii == 0                            
                    % 1. value-to-value
                    K(mask1, :) = sf2*(1 + sqrt3_r_inv_l).*exp(-sqrt3_r_inv_l);
                else
                    % 2. deri-to-value
                    K(mask1, :) = (-3*sf2*inv_ell2)*bsxfun(@minus, x(mask1, ii+1), z(:, ii)').*exp(-sqrt3_r_inv_l);
                end
            end
        end
    end
end


% derivative w.r.t log sf
if deri_wrt_param == 2
    K = 2*K;
end