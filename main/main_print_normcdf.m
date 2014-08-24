clc
clear all
close all

%% path
addpath('../util/');

X = [-8:0.1:8];
P = normcdf(X,0,1);

print_matrix_for_reference(X, 'X');
print_matrix_for_reference(P, 'P');

