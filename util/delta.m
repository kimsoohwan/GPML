function C = delta(A, B, i)
% pairwise differences w.r.t coord
%
% A: nxd
% A = [a1;  = [a11, ..., a1i, ..., a1d;
%      a2;     a21, ..., a2i, ..., a2d;
%      ...               ...
%      an];    an1, ..., ani, ..., and];
%
% B: mxd
% B = [b1;  = [b11, ..., b1i, ..., b1d;
%      b2;     b21, ..., b2i, ..., b2d;
%      ...               ...
%      bm];    bm1, ..., bmi, ..., bmd];
%
% i: coordinate
%
% C: nxm
% C = [a1i-b1i, a1i-b2i, ... a1i-bmi;
%      a2i-b1i, a2i-b2i, ... a2i-bmi;
%                ...
%      ani-b1i, ani-b2i, ... ani-bmi];
%
% Please refer to GPML/util/sq_dist.m

C = bsxfun(@minus, A(:, i), B(:, i)');