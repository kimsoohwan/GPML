function K = covMaterniso3DerObs(hyp, x, z, i)

if nargin<1, K  = '2'; return;	end     % report number of parameters
if nargin<3, z  = [];           end     % make sure, z exists
if nargin<4, i  = 0;         	end     % derivative index
dg = strcmp(z,'diag') && numel(z)>0;    % determine mode

% x, xd
if dg
    xd = [];
else
    % assume derivative observations are repeated d times
    d = size(x, 2) - 1; % first column = index
    xd_old = [];
    for dd = 1:d
        mask = x(:, 1) == d;
        xd = x(mask, 2:end);
        if dd > 1            
            assert(~any(xd_old(:) - xd(:)));
        end
        xd_old = xd;
    end
    
    xd = x(x(:, 1) == 1, 2:end);
    x  = x(x(:, 1) == 0, 2:end);
end

% call FDI template
f_handles.FF = {@covMaterniso, 3};
f_handles.FD = {@covMaterniso3FD};
f_handles.DD = {@covMaterniso3DD};
K = covDerObs(f_handles, hyp, x, xd, z, i);