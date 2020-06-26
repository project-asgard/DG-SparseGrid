clear
% clc
% close all

format short e
addpath(genpath(pwd))

tInt = 0;
tEnd = 1;

xInt = 0;
xEnd = 1;

yInt = 0;
yEnd = 1;

x_bcL = 0; x_bcR = 1;
t_bcL = 0; t_bcR = 1;


Lev = 4;
Deg = 2;
num_plot = 2;

DoFs = 2^Lev*Deg;


xmax = xEnd-xInt;
dx = xmax/2^Lev;

tmax = tEnd-tInt;
dt = tmax/2^Lev;

% Plotting Data
[quad_x,quad_w]=lgwt(num_plot,-1,1);
ww = repmat(quad_w,2^Lev,1)*dx/2;
% [x_node,Meval] = PlotDGData(Lev,Deg,xInt,xEnd,num_plot);

b = [1,1,1];


% Term 1: d/dx[FunCoef*f]
Mat_Term1_t = MatrixGrad(Lev,Deg,tInt,tEnd,1,@(x)1, @(x)0,t_bcL,t_bcR);
% Mat_Term1_t(1,1) = 1;
Mat_Term1_x = speye(DoFs,DoFs);
Mat_Term1_y = speye(DoFs,DoFs);
Mat_Term1 = kron(kron(Mat_Term1_t,Mat_Term1_x),Mat_Term1_y);

Mat_Term2_x = MatrixGrad(Lev,Deg,xInt,xEnd,1,@(x)1,@(x)0,x_bcL,x_bcR);
% Mat_Term2_x(1,1) = 1;
Mat_Term2_t = speye(DoFs,DoFs);
Mat_Term2_y = speye(DoFs,DoFs);
Mat_Term2 = kron(kron(Mat_Term2_t,Mat_Term2_x),Mat_Term2_y);

Mat_Term3_x = speye(DoFs,DoFs);
% Mat_Term2_x(1,1) = 1;
Mat_Term3_t = speye(DoFs,DoFs);
Mat_Term3_y = MatrixGrad(Lev,Deg,xInt,xEnd,1,@(x)1,@(x)0,x_bcL,x_bcR);
Mat_Term3 = kron(kron(Mat_Term3_t,Mat_Term3_x),Mat_Term3_y);

Mat = b(1)*Mat_Term1+b(2)*Mat_Term2+b(3)*Mat_Term3;
% Index_bc = [[1:DoFs],[DoFs+1:DoFs:DoFs^2]];
% for T = 1:size(Index_bc,2)
%     i = Index_bc(T);
%     Mat(i,:) = 0;
%     Mat(i,i) = 1;
% end

condest(Mat)
dof = size(Mat,1);
sol = zeros(dof,1);
rhs = zeros(dof,1);

% Exa0 = @(x,t)sin(pi*x);
% F0 = ComputRHS2D(Lev,Deg,[tInt xInt ],[tEnd xEnd ],Exa0,0);
f0 = ComputRHS(Lev,Deg,tInt,tEnd,0);
% F0 = kron(f0,[ones(Deg,1);zeros(DoFs-Deg,1)])/sqrt(dx);

% F0 = kron(kron(f0*b(1),b(2)*[legendre(-1,Deg)/sqrt(dx),zeros(1,DoFs-Deg)]'),b(3)*[legendre(-1,Deg)/sqrt(dx),zeros(1,DoFs-Deg)]');
% bc = ComputeBC(Lev,Deg,xInt,xEnd,@(x,t)sin(pi*x),0,0,0,1);

% F0 = kron(b(1)*[legendre(-1,Deg)/sqrt(dx),zeros(1,DoFs-Deg)]',kron(f0*b(1),f0*b(2)));
F0 = kron(kron(f0*b(1),f0*b(2)),b(3)*[legendre(-1,Deg)/sqrt(dx),zeros(1,DoFs-Deg)]');

[x_node,Meval] = PlotDGData(Lev,Deg,xInt,xEnd,num_plot);
[y_node,Meval] = PlotDGData(Lev,Deg,yInt,yEnd,num_plot);
[t_node,teval] = PlotDGData(Lev,Deg,tInt,tEnd,num_plot);
[x_3D_plot,y_3D_plot,z_3D_plot] = meshgrid(t_node,x_node,y_node);
MM = kron(kron(teval,Meval),Meval);
nz = size(x_3D_plot,1);

val_plot = reshape(MM*(F0),nz,nz,nz);
vtkwrite(['Initial3DLev_',num2str(Lev),'Deg_',num2str(Deg),'.vtk'],'structured_grid',x_3D_plot,y_3D_plot,z_3D_plot,...
    'scalars','sol',val_plot);
% surf(x_2D_plot,y_2D_plot,val_plot,val_plot)
% shading interp

sol = Mat\F0;
% figure
val_plot = reshape(MM*(sol),nz,nz,nz);
% % surf(x_2D_plot,y_2D_plot,val_plot,val_plot)
% scatter3(x_3D_plot,y_3D_plot,z_3D_plot,val_plot)
% shading interp

vtkwrite(['LinearConvect3DLev_',num2str(Lev),'Deg_',num2str(Deg),'.vtk'],'structured_grid',x_3D_plot,y_3D_plot,z_3D_plot,...
    'scalars','sol',val_plot);