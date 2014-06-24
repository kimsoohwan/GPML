function TEST_EQ(A, B, test_name, eps)

% epsilon
if nargin < 4
    eps = 10^(-5);
end

% test NaN
if any(isnan(A(:)))
    fprintf(2, '\t\tWarning: A has NaN\n');
end
if any(isnan(B(:)))
    fprintf(2, '\t\tWarning: B has NaN\n');
end

% test equality
if any(abs(A(:) - B(:)) >= eps)
    fprintf(2, ['\t\tFailure: ', test_name, '\n']);
else
    fprintf(1, ['\t\tSuccess: ', test_name, '\n']);
end