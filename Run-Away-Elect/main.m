% This code tests the run-away electron problem

clear
close all
format short e

addpath(genpath(pwd))

Lev = 2;
Deg = 3;
Lmax = 1;
pde = RunAwayElect1;

% Mat1D = Matrix(Lev,Deg,Lmax);
% Rhs1D = Rhs(Deg,Lev,Lmax,pde);

p_1 = legendre(-1,Deg);
p_2 = legendre( 1,Deg);

quad_num = 10;

[quad_x,quad_w] = lgwt(quad_num,-1,1);
p_val = legendre(quad_x,Deg);
Dp_val = dlegendre(quad_x,Deg);
Dp_1 = dlegendre(-1,Deg);
Dp_2 = dlegendre(1,Deg);


%---------------------------
% Define dofs
%---------------------------
nx = 2^(Lev);
hx = Lmax/nx;
dofs = Deg*nx;
