function print_matrix_for_reference(K, strK)

disp([strK, ': ', num2str(size(K, 1)), ' x ' , num2str(size(K, 2))]);
for row = 1:size(K, 1)
    for col = 1:size(K, 2)
        if K(row, col) >= 0
            fprintf(' %.15ff, ', K(row, col));
        else
            fprintf('%.15ff, ', K(row, col));
        end
    end
    fprintf('\n');
end
fprintf('\n');