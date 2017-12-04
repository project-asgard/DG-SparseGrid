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
% return

%% 2D sparse-grid
tic
run Compute_2D_sparsegrid_v3.m
toc
['Compute_2D_sparsegrid.m is done']

tic
run Compute_4D_sparsegrid.m
toc
['Compute_2D_sparsegrid.m is done']


function I=Index1D(mx,ix,n)
% index for 1D
I=(2^mx+ix);
end



