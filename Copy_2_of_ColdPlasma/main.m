%--------------------------------------------------------------
% This code solves the Maxwell's equations by DG-SG methods
% E=[E1(x,y,z),E2(x,y,z),E3(x,y,z)]^t
%--------------------------------------------------------------
clc
clear
close all

format short e

addpath(genpath(pwd))

%% Step 1. Setting Parameters
Lev = 2;
Deg = 2;
Dim = 3;

Lmax = 1;
pde = Maxwell1;

%*************************************************
%% Step 1.1. Set Up Matrices for Multi-wavelet
% Input:  Deg and Lev
% Output: Convert Matrix FMWT_COMP
%*************************************************
FMWT_COMP_x = OperatorTwoScale(Deg,2^Lev);

%% Coefficient Matrix for time-independent matrix
[GradMat,MaxMat] = Matrix_TI(Lev,Deg,Lmax,FMWT_COMP_x);

return
%1D connectivity
Con1D=Connect1D(Lev);
%% Hash Table
% Sparse Grids Simulation
[Hash,InvHash]=HashTable(Lev,Dim);

%% Global Matrix for 3D Maxwell Eq.
% global rhs, E0, and B0 vectors
[b_s,E_s,B_s]=GlobalRHS(Deg,F_1D,E_1D,B_1D,InvHash);

%% Maxwell Solver
[Eh,Bh] = MaxwellSolver(Lev,Deg,Hash,InvHash,Con1D,GradX,pde.eps,pde.mu,pde.w,dt,MaxT,b_s,E_s*cos(0),B_s*0);
sol_n=[Eh;Bh];

% [Eh,Bh] = MaxwellSolver2(Lev,Deg,Hash,InvHash,Con1D,GradX,pde.eps,pde.mu,pde.w,dt,MaxT,b_s,E_s*cos(0),B_s*0);
% sol_n=[Eh;Bh];

%% Error Estimate
time=dt*MaxT;
u_s=[E_s*cos(pde.w*time);B_s*sin(pde.w*time)];

full([Deg Lev max(abs(sol_n-u_s)) norm(sol_n-u_s)])