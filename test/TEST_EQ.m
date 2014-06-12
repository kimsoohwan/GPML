function TEST_EQ(A, B, test_name, eps)

% epsilon
if nargin < 4
    eps = 10^(-5);
end

% test equality
if any(abs(A - B) >= eps)
    fprintf(1, ['\t\tFailure: ', test_name, '\n']);
else
    fprintf(1, ['\t\tSuccess: ', test_name, '\n']);
end