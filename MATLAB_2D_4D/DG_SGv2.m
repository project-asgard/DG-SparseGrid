%%Test Multiwavelet Sparse Grids for DG scheme%%
% 2D Elliptic Problem
% 4D Elliptic Problem

close all
clear all
clc
format short e


%------------------------------------------------
% Set Parameters
%------------------------------------------------
Np = 10;
k=5;

mkdir('Data')
addpath(genpath(pwd))

% n: the level for finest mesh
% l: places for the functions, l=0,...,2^n-1
% k: the degrees of polynomials (k-1)

n = Np;h = 2^(-n-1);
%% Compute 1D solution
[Stiff_1D,b,Meval,coef_MW] = LaplacianMatrix2(n,k);
M_mass=speye(size(Stiff_1D,1),size(Stiff_1D,1));

%% 2D sparse-grid
tic
run Compute_2D_sparsegridv2.m
% run Compute_2D_sparsegrid.m
toc
['Compute_2D_sparsegridv2.m is done']

return
tic
run Compute_4D_sparsegridv2.m
% run Compute_4D_sparsegrid.m
toc
['Compute_4D_sparsegridv2.m is done']


%  function I=Index1D(mx,ix,n)
%  % index for 1D
%  I=(2^mx+ix);
%  end
