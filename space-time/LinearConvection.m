% Implementation for solve 1D space-time 
% convection equation

clear
clc
close all

format short e
addpath(genpath(pwd))

tInt = 0;
tEnd = 1;

xInt = 0;
xEnd = 1;

x_bcL = 0; x_bcR = 1;
t_bcL = 0; t_bcR = 1;


Lev = 5;
Deg = 4;
num_plot = 2;

DoFs = 2^Lev*Deg;


xmax = xEnd-xInt;
dx = xmax/2^Lev;

tmax = tEnd-tInt;
dt = tmax/2^Lev;

% Define the mesh
Dim = 2;
[database,Inv] = HashTable(Lev,Lev,Deg,Dim);
DoF_SG = size(Inv.x1,2);
Mat_SG = sparse(DoF_SG,DoF_SG);

% Plotting Data
[quad_x,quad_w]=lgwt(num_plot,-1,1);
ww = repmat(quad_w,2^Lev,1)*dx/2;


b = [1,1];

FMWT_COMP = OperatorTwoScale(Deg,2^Lev);

% Term 1: d/dt[f]
% DG basis
Mat_Term1_t = MatrixGrad(Lev,Deg,tInt,tEnd,1,@(x)1, @(x)0,t_bcL,t_bcR);
% Wavelet basis
Mat_Term1_t = FMWT_COMP*Mat_Term1_t*FMWT_COMP';
Mat_Term1_x = speye(DoFs,DoFs);
% Assemble Term1 to Sparse Grids Matrix

% for i = 1:DoF_SG
%     for j = 1:DoF_SG
%         I1 = Inv.x1(i);
%         I2 = Inv.x2(i);
%         Ind1 = Inv.x1(j);
%         Ind2 = Inv.x2(j);
%         val = Mat_Term1_t(I1,Ind1)*Mat_Term1_t(I2,Ind2);
%         Mat_SG(i,j) = Mat_SG(i,j)+val;
%     end
% end


% Full grids
Mat_Term1 = kron(Mat_Term1_t,Mat_Term1_x);

% Term 2: d/dx[FunCoef*f]
Mat_Term2_x = MatrixGrad(Lev,Deg,xInt,xEnd,1,@(x)1,@(x)0,x_bcL,x_bcR);
Mat_Term2_t = speye(DoFs,DoFs);
Mat_Term2 = kron(Mat_Term2_t,Mat_Term2_x);

Mat = b(1)*Mat_Term1+b(2)*Mat_Term2;


dof = size(Mat,1);
sol = zeros(dof,1);
rhs = zeros(dof,1);

% Exa0 = @(x,t)sin(pi*x);
% F0 = ComputRHS2D(Lev,Deg,[tInt xInt ],[tEnd xEnd ],Exa0,0);
f0 = ComputRHS(Lev,Deg,tInt,tEnd,0);
F0 = kron(f0*b(1),b(2)*[legendre(-1,Deg)/sqrt(dx),zeros(1,DoFs-Deg)]');


[x_node,Meval] = PlotDGData(Lev,Deg,xInt,xEnd,num_plot);
[t_node,teval] = PlotDGData(Lev,Deg,tInt,tEnd,num_plot);
[x_2D_plot,y_2D_plot] = meshgrid(t_node,x_node);
MM = kron(teval,Meval);
nz = size(x_2D_plot,1);

val_plot = reshape(MM*(F0),nz,nz);
% surf(x_2D_plot,y_2D_plot,val_plot,val_plot)
% shading interp

sol = Mat\F0;
figure
val_plot = reshape(MM*(sol),nz,nz);
surf(x_2D_plot,y_2D_plot,val_plot,val_plot)
shading interp
colormap jet