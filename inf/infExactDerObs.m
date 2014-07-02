function [post nlZ dnlZ] = infExactDerObs(hyp, mean, cov, lik, x, y)

% Exact inference for a GP with Gaussian likelihood. Compute a parametrization
% of the posterior, the negative log marginal likelihood and its derivatives
% w.r.t. the hyperparameters. See also "help infMethods".
%
% Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch, 2011-02-18
%
% See also INFMETHODS.M.

likstr = lik; if ~ischar(lik), likstr = func2str(lik); end 
if ~strcmp(likstr,'likGaussDerObs')         % NOTE: no explicit call to likGauss
  error('Exact inference only possible with Gaussian likelihood dealing with derivative observations');
end

% number of function and derivative training data
nd = sum(x(:, 1) ~= 0);     % derivative training data
n = size(x, 1) - nd;        % function derivative training data
d = size(x, 2) - 1;         % number of dimensions
nn = n + nd*d;

K = feval(cov{:}, hyp.cov, x);                      % evaluate covariance matrix
m = feval(mean{:}, hyp.mean, x);                          % evaluate mean vector

sn2 = exp(2*hyp.lik(1));                            % noise variance of likGauss
snd2 = exp(2*hyp.lik(2));                           % noise variance of likGauss
D = [ sn2*ones(n, 1);
      snd2*ones(nd*d, 1)];
L = chol(K+diag(D));               % Cholesky factor of covariance with noise
alpha = solve_chol(L,y-m);

post.alpha = alpha;                            % return the posterior parameters
post.sW = ones(nn,1);                  % sqrt of noise precision vector
post.L  = L;                                        % L = chol(eye(n)+sW*sW'.*K)

if nargout>1                               % do we want the marginal likelihood?
    % p(y) = N(m, Kn) = (2pi)^(-n/2) * |Kn|^(-1/2) * exp[(-1/2) * (y-m)' * inv(Kn) * (y-m)]
    % nlZ = (1/2) * (y-m)' * inv(Kn) * (y-m) + (1/2) * log |Kn|		+ (n/2) * log(2pi)
    %     = (1/2) * (y-m)' * alpha           + (1/2) * log |L*L'|		+ (n/2) * log(2pi)
    %     = (1/2) * (y-m)' * alpha           + (1/2) * log |L|*|L'|	+ (n/2) * log(2pi)
    %     = (1/2) * (y-m)' * alpha           + log |L||					+ (n/2) * log(2pi)
    %     = (1/2) * (y-m)' * alpha           + tr[log (L)]				+ (n/2) * log(2pi)
	nlZ = (y-m)'*alpha/2 + sum(log(diag(L))) + n*log(2*pi)/2;  % -log marg lik
    
    if nargout>2                                         % do we want derivatives?
        dnlZ = hyp;                                 % allocate space for derivatives
        
        % (1) w.r.t the mean parameters
        % nlZ = (1/2) * (y-m)' * inv(Kn) * (y-m)
        %       = - m' * inv(Kn) * y + (1/2) m' * inv(Kn) * m
        % nlZ_i = - m_i' * inv(Kn) * y + m_i' * inv(Kn) * m
        %       = - m_i' * inv(Kn) (y - m)
        %       = - m_i' * alpha
        for i = 1:numel(hyp.mean), 
            dnlZ.mean(i) = -feval(mean{:}, hyp.mean, x, i)'*alpha;
        end
        
        % (2) w.r.t the cov parameters
        % nlZ = (1/2) * (y-m)' * inv(Kn) * (y-m) + (1/2) * log |Kn|
        % nlZ_j = (-1/2) * (y-m)' * inv(Kn) * K_j * inv(Kn) * (y-m)	+ (1/2) * tr[inv(Kn) * K_j]
        %       = (-1/2) * alpha' * K_j * alpha						+ (1/2) * tr[inv(Kn) * K_j]
        %       = (-1/2) * tr[(alpha' * alpha) * K_j]				+ (1/2) * tr[inv(Kn) * K_j]
        %       = (1/2) tr[(inv(Kn) - alpha*alpha') * K_j]
        %       = (1/2) tr[Q * K_j]
        % Q = inv(Kn) - alpha*alpha'
        Q = solve_chol(L,eye(nn)) - alpha*alpha';    % precompute for convenience
        for i = 1:numel(hyp.cov)
            dnlZ.cov(i) = sum(sum(Q.*feval(cov{:}, hyp.cov, x, [], i)))/2;
        end
 
        % (3) w.r.t the cov parameters
        % nlZ = (1/2) * (y-m)' * inv(K + D) * (y-m) + (1/2) * log |K + D|
        % nlZ_j = (-1/2) * (y-m)' * inv(Kn) * D_j * inv(Kn) * (y-m)	+ (1/2) * tr[inv(Kn) * D_j]
        %       = (-1/2) * alpha' * D_j * alpha						+ (1/2) * tr[inv(Kn) * D_j]
        %       = (-1/2) * tr[(alpha' * alpha) * D_j]				+ (1/2) * tr[inv(Kn) * D_j]
        %       = (1/2) tr[(inv(Kn) - alpha*alpha') * D_j]
        %       = (1/2) tr[Q * D_j]       
        dnlZ.lik(1) = sn2*trace(Q(1:n, 1:n));
        dnlZ.lik(2) = snd2*trace(Q(n+1:end, n+1:end));
    end
end
