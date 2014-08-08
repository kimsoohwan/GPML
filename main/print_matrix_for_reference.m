function print_matrix_for_reference(K, strK, num_col)
% number of columns
if nargin < 3
    num_col = size(K, 2);
end

disp([strK, ': ', num2str(size(K, 1)), ' x ' , num2str(size(K, 2))]);
col_index = 0;
for row = 1:size(K, 1)
    for col = 1:size(K, 2)
        if K(row, col) >= 0
            fprintf(' %.15ff', K(row, col));
        else
            fprintf('%.15ff', K(row, col));
        end
        col_index = col_index + 1;
        if col_index == num_col
            if row ~= size(K, 1)
                fprintf(',\n');
                col_index = 0;
            else
                fprintf(';\n')
            end
        else
            fprintf(', ');
        end
    end
    %fprintf('\n');
end
fprintf('\n');