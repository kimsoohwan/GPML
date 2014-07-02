function A = meanZeroDerObs(hyp, x, i)

% Zero mean function. The mean function does not have any parameters.
%
% m(x) = 0

if nargin<2, A = '0'; return; end             % report number of hyperparameters 

% number of function and derivative training data
nd = sum(x(:, 1) ~= 0);     % derivative training data
n = size(x, 1) - nd;        % function derivative training data
d = size(x, 2) - 1;         % number of dimensions
nn = n + nd*d;

% value
A = zeros(nn, 1);           % derivative and mean
