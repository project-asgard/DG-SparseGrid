%===================================================
% This code solves Maxwell's equation by
% DG-SG grids methods
%===================================================
%clc
clear
close all

format short e

addpath(genpath(pwd))

%% Step 1. Setting Parameters
Lev = 4;
Deg = 2;


Lmax = 1;
pde = Maxwell1;
% cfl=4000/2/2/2/2/2/2/2;
dt = 8*10^(-4);
%dt=1/80;
MaxT =ceil(0.5/dt);%10;



%*************************************************
%% Step 1.1. Set Up Matrices for Multi-wavelet
% Input:  Deg and Lev
% Output: Convert Matrix FMWT_COMP
%*************************************************
FMWT_COMP_x = OperatorTwoScale(Deg,2^Lev);

%*************************************************
%% Step 1.2. Set Up Initial Conditions
% F_1D: rhs coefficients
% E_1D: Initial Condition for E
% B_1D: Initial condition for B
% Input:  Deg, Lev,Lmax, pde,FMWT_COMP_x
% Output: F_1D,E_1D,B_1D
%*************************************************
[F_1D,E_1D,B_1D] = Intial_Con(Lev,Deg,Lmax,pde,FMWT_COMP_x);

%*************************************************
%% Step 2. Hash Table and 1D Connectivity
%*************************************************
[Hash,InvHash]=HashTable(Lev);
%1D connectivity
Con1D=Connect1D(Lev);

%*************************************************
%% Step 3. Coefficient Matrix for Time-independent Matrix
%*************************************************
GradX = Matrix_TI(Lev,Deg,Lmax,FMWT_COMP_x);

%*************************************************
%% Step 4.3D Maxwell Solver
%*************************************************
% global rhs, E0, and B0 vectors
% b_s is the RHS vector
% E_s and B_s are used for error estimate

%% Maxwell Solver
[Bh,E1h,E2h] = MaxwellSolver1(Lev,Deg,Hash,InvHash,Con1D,GradX,pde.w,dt,MaxT,...
    F_1D,E_1D.E1,E_1D.E2,B_1D.B);
sol_n=[Bh;E1h;E2h];

%% Error Estimate
time=dt*MaxT;
u_s=[B_1D.B*exp(-time);E_1D.E1*exp(-time);E_1D.E2*exp(-time)];

Bs=[B_1D.B*exp(-time)];E1s=[E_1D.E1*exp(-time)];E2s=[E_1D.E2*exp(-time)];

% u_s=[B_1D.B;E_1D.E1;E_1D.E2];
% 
% Bs=[B_1D.B];E1s=[E_1D.E1];E2s=[E_1D.E2];
full([Deg Lev dt max(abs(sol_n-u_s)) norm(sol_n-u_s)])

% subplot(1,2,1)
% plot(sol_n)
% subplot(1,2,2)
% plot(u_s)
% 
% [M,N]=matrix_plot(Lev,Deg,Lmax,FMWT_COMP_x);
% figure;
% subplot(1,3,1)
% plot(N,M*E1h,'r--',N,M*E1s,'b--')
% legend('Numerical Sol','Real Solution')
% title('E1')
% subplot(1,3,2)
% plot(N,M*E2h,'r--',N,M*E2s,'b--')
% legend('Numerical Sol','Real Solution')
% title('E2')
% subplot(1,3,3)
% plot(N,M*Bh,'r--',N,M*Bs,'b--')
% legend('Numerical Sol','Real Solution')
% title('B')

