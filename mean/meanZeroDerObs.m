function A = meanZeroDerObs(hyp, x, i)

% Zero mean function. The mean function does not have any parameters.
%
% m(x) = 0

if nargin<2, A = '0'; return; end             % report number of hyperparameters 

% number of function and derivative training data
[nn, d] = size(x);
d = d - 1; % first column = index

% value
A = zeros(nn, 1);           % derivative and mean
