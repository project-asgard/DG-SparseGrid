%%Test Multiwavelet for DG scheme%%
% 2D Elliptic Problem

close all
clear all
clc
format short e
%---------------------------------------------
% ?????
% 1. needs more efficient strategy for MD-Index
%---------------------------------------------

%------------------------------------------------
% Set Parameters
%------------------------------------------------
Np = 5;
k=3;
scheme='sparse';
mkdir('Data')

% n: the level for finest mesh
% l: places for the functions, l=0,...,2^n-1
% k: the degrees of polynomials (k-1)

n = Np;h = 2^(-n-1);
%% Compute 1D solution
run Comp_1D_DG.m
['Comp_1D_DG.m is done']

return

%% 2D full-grid
run Compute_2D_fullgrid.m
['Compute_2D_fullgrid.m is done']

%% 2D sparse-grid
run Compute_2D_sparsegrid.m
['Compute_2D_sparsegrid.m is done']

%% comparing sparsity
run Sparsity.m
['Sparsity.m is done']

run PlotSolution.m
['PlotSolution.m is done']

% Comparing Error and Condition Number for Sparse Grid and Full Grid
disp(' ')
disp('Max Error of 2D Computing     Condition Number of 2D Computing ')
disp('Sparse Grid       Full Grid           Sparse Grid       Full Grid')
formatSpec = '%6.4e      %6.4e      %6.4e      %6.4e \n';
fprintf(formatSpec, Lmax_2D_s,Lmax_2D,cond_sparse,cond_full)


function I=Index1D(mx,ix,n)
% index for 1D
I=(2^mx+ix);
end



