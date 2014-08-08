clc
clear all
close all

%% path
addpath('../util/');


%% setting
min_pt.x = 0.5377;
min_pt.y = 1.8339;
min_pt.z = -2.2588;
n = 3; % number of cells per axis
cell_size = 0.1;
half_cell_size = cell_size/2;

%% test data
xs = zeros(n*n*n, 3);
idx = 1;
for i = 1:n
    x = min_pt.x + half_cell_size + (i-1)*cell_size;    
    for j = 1:n
        y = min_pt.y + half_cell_size + (j-1)*cell_size;
        for k = 1:n
            z = min_pt.z+ half_cell_size + (k-1)*cell_size;
            xs(idx, 1) = x;
            xs(idx, 2) = y;
            xs(idx, 3) = z;
            idx = idx+1;
        end
    end
end

print_matrix_for_reference(xs, 'xs');