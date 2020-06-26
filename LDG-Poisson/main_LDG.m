% main code for LDG Poisson equation
format short e
addpath(genpath(pwd))

Lev = 2;
Deg = 1;
Lmax = 1;

%--Quadrature
quad_num=10;
%---------------

% compute the trace values
p_1 = legendre(-1,Deg);
p_2 = legendre(1,Deg);

[quad_x,quad_w]=lgwt(quad_num,-1,1);
p_val = legendre(quad_x,Deg);
Dp_val = dlegendre(quad_x,Deg);

%---------------------------
% Jacobi of variable x and v
% Define Matrices
%---------------------------
nx=2^(Lev);hx=Lmax/nx;
Jacobi_x=hx;
dof_1D_x=Deg*nx;
GradX=sparse(dof_1D_x,dof_1D_x);
DeltaX=sparse(2*dof_1D_x,2*dof_1D_x);

FMWT_COMP_x = OperatorTwoScale(Deg,2^Lev);

%***************************
% Term for \int_K(u_h)(v')dK
%***************************
val = 1/hx*[Dp_val'*(quad_w.*p_val)];
Ac = repmat({val},nx,1);
GradX = blkdiag(Ac{:});

%****************************************
% Term for \int_{\partial K} {{u_h}}vds
% Numerical Flux is taken as central flux
%****************************************
Amd  = -p_1'*p_1/2/hx+p_2'*p_2/2/hx;
Asub = -p_1'*p_2/2/hx;
Asup = zeros(Deg,Deg)+p_2'*p_1/2/hx;
GradXFluxC = -blktridiag([Amd],[Asub],[Asup],nx);
% Adding Periodic Boundary Conditions
IndexEnd = dof_1D_x-Deg+1:dof_1D_x;
IndexSta = 1:Deg;
Iu = repmat(IndexEnd,Deg,1);
Iv = repmat(IndexSta,Deg,1);
GradXFluxC = GradXFluxC...
    +sparse([Iv',Iu'],[Iu,Iv],-[Asub,Asup],dof_1D_x,dof_1D_x);

%**************************************
% Term for \int_{\partial K} [[u_h]]vds
% Numerical Flux is taken as jump
%**************************************
Amd  = -p_1'*p_1/hx+p_2'*(-p_2)/hx;
Asub = -p_1'*(-p_2)/hx;
Asup = zeros(Deg,Deg)+p_2'*p_1/hx;
GradXFluxJ = -blktridiag([Amd],[Asub],[Asup],nx);

% Adding Periodic Boundary Conditions
IndexEnd = dof_1D_x-Deg+1:dof_1D_x;
IndexSta = 1:Deg;
Iu = repmat(IndexEnd,Deg,1);
Iv = repmat(IndexSta,Deg,1);
GradXFluxJ = GradXFluxJ...
    +sparse([Iv',Iu'],[Iu,Iv],-[Asub,Asup],dof_1D_x,dof_1D_x);


% PGradX denote the flux is taken as \hat{u}=u+
% NGradX denote the flux is taken as \hat{u}=u-
% GradX denote the flux is taken as \hat{u}={u}
PGradX = GradX+GradXFluxC+GradXFluxJ/2;
NGradX = GradX+GradXFluxC-GradXFluxJ/2;
 GradX = GradX+GradXFluxC;
 
 
 

A_Poisson = [speye(dof_1D_x) PGradX;...
    NGradX sparse(dof_1D_x,dof_1D_x)];


Mat_S = NGradX*PGradX;
Mat_S(1,:)=0;
Mat_S(1,1)=1;
% Mat_S(end,:)=0;
% Mat_S(end,end)=1;
A_Poisson(dof_1D_x+1,:)=0;
A_Poisson(dof_1D_x+1,dof_1D_x+1)=1;

DeltaX = blkdiag(FMWT_COMP_x,FMWT_COMP_x)*...
                A_Poisson*...
         blkdiag(FMWT_COMP_x',FMWT_COMP_x');
FMat_S =   FMWT_COMP_x*Mat_S*FMWT_COMP_x;   
     
% A_Poisson(end,:)=0;A_Poisson(end,end)=1;
[condest(Mat_S) condest(FMat_S) condest(A_Poisson) condest(DeltaX)]
[size(Mat_S) size(A_Poisson)]



figure;
subplot(1,2,1);spy(A_Poisson)
subplot(1,2,2);spy(Mat_S)
% (Mat_S)\ones(dof_1D_x,1)
return
%---------------------------
% Jacobi of variable x and v
%---------------------------
nx=2^(Lev);hx=Lmax/nx;
dof_1D_x=k*nx;

b_poisson=sparse(2*dof_1D_x,1);


%-----------------------------------------------------
% Handling B.C. for enforcing phi(0)=phi(Lmax)=0
%-----------------------------------------------------
% Phi_v(0,0) corresponding to the

Index=[1:k]'+[0:k^2:(2^(Lev_x)-1)*k^2];
tmp_b=[sqrt(Lmax); zeros(dof_1D_x-1,1)]...
      -sqrt(2*Vmax)*f(Index(:));%f(1:dof_1D_x);


b_poisson(dof_1D_x+1:end)=(FMWT_COMP_x(:,2:end-1)*FMWT_COMP_x(:,2:end-1)')*tmp_b;


x_poisson=DeltaX\b_poisson;

EE=x_poisson(1:dof_1D_x);