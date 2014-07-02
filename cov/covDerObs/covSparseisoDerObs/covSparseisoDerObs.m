function K = covSparseisoDerObs(hyp, x, z, i)

if nargin<1, K  = '2'; return;	end     % report number of parameters
if nargin<3, z  = [];           end     % make sure, z exists
if nargin<4, i  = 0;         	end     % derivative index
dg = strcmp(z,'diag') && numel(z)>0;    % determine mode

% x, xd
if dg
    xd = [];
else
    der_mask = x(:, 1) ~= 0;
    xd = x( der_mask, 2:end);
    x  = x(~der_mask, 2:end);
end

% call FDI template
f_handles.FF = {@covSparseiso};
f_handles.FD = {@covSparseisoFD};
f_handles.DD = {@covSparseisoDD};
K = covDerObs(f_handles, hyp, x, xd, z, i);