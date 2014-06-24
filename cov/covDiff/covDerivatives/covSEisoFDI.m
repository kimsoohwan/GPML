function K = covSEisoFDI(n2, hyp, x, z, ii)

if nargin<2, K = '2'; return; end                  % report number of parameters
if nargin<4, z = []; end                                   % make sure, z exists
if nargin<5, ii = 0; end                                      % derivative index

% call FDI template
hCovFF = {@covSEiso};
hCovFD = {@covSEisoFD};
hCovDD = {@covSEisoDD};
K = covFDI(n2, hyp, x, z, ii, hCovFF, hCovFD, hCovDD);