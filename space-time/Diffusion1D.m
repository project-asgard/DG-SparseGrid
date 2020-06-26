% clear
% clc
close all

format short e
addpath(genpath(pwd))

tInt = 0;
tEnd = 1;

xInt = 0;
xEnd = 1;

x_bcL = 0; x_bcR = 0;
t_bcL = 0; t_bcR = 1;


Lev = 6;
Deg = 6;
num_plot = 2;

DoFs = 2^Lev*Deg;


xmax = xEnd-xInt;
dx = xmax/2^Lev;

tmax = tEnd-tInt;
dt = tmax/2^Lev;

% Plotting Data
[quad_x,quad_w]=lgwt(num_plot,-1,1);
ww = repmat(quad_w,2^Lev,1)*dx/2;


% b = [pi^2,-1];

b = [pi^2,-1];

% Term 1: d/dx[FunCoef*f]
Mat_Term1_t = MatrixGrad(Lev,Deg,tInt,tEnd,1,@(x)1, @(x)0,t_bcL,t_bcR);
Mat_Term1_x = speye(DoFs,DoFs);
Mat_Term1 = kron(Mat_Term1_t,Mat_Term1_x);

tmp1 = MatrixGrad(Lev,Deg,xInt,xEnd,0,@(x)1,@(x)0,1,1);
tmp2 = MatrixGrad(Lev,Deg,xInt,xEnd,0,@(x)1,@(x)0,0,0);

% tmp1 = MatrixGrad(Lev,Deg,xInt,xEnd,1,@(x)1-x.^2,@(x)0,1,1);
% tmp2 = MatrixGrad(Lev,Deg,xInt,xEnd,-1,@(x)1,@(x)0,0,0);
Mat_Term2_x = tmp1*tmp2;
Mat_Term2_t = speye(DoFs,DoFs);
Mat_Term2 = kron(Mat_Term2_t,Mat_Term2_x);

Mat = b(1)*Mat_Term1+b(2)*Mat_Term2;


dof = size(Mat,1);
sol = zeros(dof,1);
rhs = zeros(dof,1);


f0 = ComputRHS(Lev,Deg,xInt,xEnd,tInt);
% f1 = ComputRHS(Lev,Deg,tInt,tEnd,tEnd);
F0 = kron(b(1)*[legendre(-1,Deg)/sqrt(dt),zeros(1,DoFs-Deg)]',-b(2)*f0);


[x_node,Meval] = PlotDGData(Lev,Deg,xInt,xEnd,num_plot);
[t_node,teval] = PlotDGData(Lev,Deg,tInt,tEnd,num_plot);
[x_2D_plot,y_2D_plot] = meshgrid(x_node,t_node);
MM = kron(teval,Meval);
nz = size(x_2D_plot,1);

% val_plot = reshape(MM*(F0),nz,nz);
% surf(x_2D_plot,y_2D_plot,val_plot,val_plot)
% shading interp


sol = Mat\F0;
figure
val_plot = reshape(MM*(sol),nz,nz);
surf(x_2D_plot,y_2D_plot,val_plot',val_plot')
shading interp
colormap jet

uu = exp(-y_2D_plot).*sin(pi*x_2D_plot);
% figure;
% surf(x_2D_plot,y_2D_plot,uu,uu)
% shading interp
% colormap jet

figure;
surf(x_2D_plot,y_2D_plot,val_plot'-uu,val_plot'-uu)
shading interp
colormap jet

err = val_plot'-uu;
err = err(:);
[Deg Lev condest(Mat) norm(err)]
