function K = covMaterniso3FDI(n2, hyp, x, z, ii)

if nargin<2, K = '2'; return; end                  % report number of parameters
if nargin<4, z = []; end                                   % make sure, z exists
if nargin<5, ii = 0; end                                      % derivative index

% call FDI template
hCovFF = {@covMaterniso, 3};
hCovFD = {@covMaterniso3FD};
hCovDD = {@covMaterniso3DD};
K = covFDI(n2, hyp, x, z, ii, hCovFF, hCovFD, hCovDD);