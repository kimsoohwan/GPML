function K = covSEisoDiffFast(hyp, x, z, i, xd)

if nargin<1, K  = '2'; return;	end     % report number of parameters
if nargin<3, z  = [];        	end  	% make sure, z exists
if nargin<4, i  = 0;        	end  	% derivative index
if nargin<5, xd = [];           end  	% derivative inputs

% call FDI template
f_handles.FF = {@covSEiso};
f_handles.FD = {@covSEisoFD};
f_handles.DD = {@covSEisoDD};
K = covDiffFast(f_handles, hyp, x, xd, z, i);