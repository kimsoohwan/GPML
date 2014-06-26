function K = covSparseisoDiff(hyp, x, z, i, xd)

if nargin<1, K  = '2'; return;	end     % report number of parameters
if nargin<3, z  = [];           end     % make sure, z exists
if nargin<4, i  = 0;         	end     % derivative index
if nargin<5, xd = [];           end  	% derivative inputs

% call FDI template
f_handles.FF = {@covSparseiso};
f_handles.FD = {@covSparseisoFD};
f_handles.DD = {@covSparseisoDD};
K = covDiff(f_handles, hyp, x, xd, z, i);