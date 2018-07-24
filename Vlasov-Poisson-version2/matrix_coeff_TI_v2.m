function [GradX]=matrix_coeff_TI_v2(Lev,k,Lmax,Vmax)
%=============================================================
% Generate time-independent coefficient matrices
% Vlasolv Solver:
%   Operators:  vMassV: int_v v*l_i(v)*l_j(v)dv
%               GradV: int_v (l_i(v))'*l_j(v)dv
%               GradX: int_x (m_i(x))'*m_j(x)dx
% Poisson Solver:
%   Operators: DeltaX: int_x (m_i(x))''*m_j(x)dx
%   This equation is solved by LDG methods
% Maxwell Solver: (Missing)
%   Operators: CurlX: int_x curl(m_i(x))*m_j(x)dx
% Input: Lev, k, dim, Lmax, Vmax
% P.S. This is the full-grid version
%=============================================================
%--Quadrature
quad_num=10;
%---------------

% compute the trace values
p_1 = legendre(-1,k);
p_2 = legendre(1,k);

[quad_x,quad_w]=lgwt(quad_num,-1,1);
p_val = legendre(quad_x,k);
Dp_val = dlegendre(quad_x,k);

%---------------------------
% Jacobi of variable x and v
% Define Matrices
%---------------------------
nx=2^(Lev);hx=Lmax/nx;
Jacobi_x=hx;
dof_1D_x=k*nx;
GradX=sparse(dof_1D_x,dof_1D_x);
DeltaX=sparse(2*dof_1D_x,2*dof_1D_x);


nv=2^(Lev);hv=2*Vmax/nv;
Jacobi_v=hv;
dof_1D_v=k*nv;
vMassV=sparse(dof_1D_v,dof_1D_v);
GradV=sparse(dof_1D_v,dof_1D_v);

%***************************
% Term for \int_K(u_h)(v')dK
%***************************
val = 1/hx*[Dp_val'*(quad_w.*p_val)];
Ac = repmat({val},nx,1);
GradX = blkdiag(Ac{:});

full(GradX)
load(['two_scale_rel_',num2str(k),'.mat'])

for j_x = 1:k
    for j_y = 1:k
        H1(j_x,j_y) = ((-1)^(j_x+j_y-2)    )*H0(j_x,j_y);
        G1(j_x,j_y) = ((-1)^(k+j_x+j_y-2))*G0(j_x,j_y);
    end
end


size(GradX)
GG=[G0 G1 zeros(2,4) zeros(2,4) zeros(2,4);...
    zeros(2,4) G0 G1 zeros(2,4) zeros(2,4);...
    zeros(2,4) zeros(2,4) G0 G1 zeros(2,4);...
    zeros(2,4) zeros(2,4) zeros(2,4) G0 G1];
size(GG)
full(GG*GradX*GG')

FMWT_COMP_x = OperatorTwoScale(k,2^Lev);


GradX  = FMWT_COMP_x*GradX*FMWT_COMP_x';
full(GradX)
% % % figure out the matrix here
% % GradX [H0 H1]*GradX*[G0 G1]'
% %   [G0 G1]*GradX*[H0 H1]' [G0 G1]*GradX*[G0 G1]'
return

%****************************************
% Term for \int_{\partial K} {{u_h}}vds
% Numerical Flux is taken as central flux
%****************************************
Amd  = -p_1'*p_1/2/hx+p_2'*p_2/2/hx;
Asub = -p_1'*p_2/2/hx;
Asup = zeros(k,k)+p_2'*p_1/2/hx;
GradXFluxC = -blktridiag([Amd],[Asub],[Asup],nx);
% Adding Periodic Boundary Conditions
IndexEnd = dof_1D_x-k+1:dof_1D_x;
IndexSta = 1:k;
Iu = repmat(IndexEnd,k,1);
Iv = repmat(IndexSta,k,1);
GradXFluxC = GradXFluxC...
    +sparse([Iv',Iu'],[Iu,Iv],-[Asub,Asup],dof_1D_x,dof_1D_x);

%**************************************
% Term for \int_{\partial K} [[u_h]]vds
% Numerical Flux is taken as jump
%**************************************
Amd  = -p_1'*p_1/hx+p_2'*(-p_2)/hx;
Asub = -p_1'*(-p_2)/hx;
Asup = zeros(k,k)+p_2'*p_1/hx;
GradXFluxJ = -blktridiag([Amd],[Asub],[Asup],nx);

% Adding Periodic Boundary Conditions
IndexEnd = dof_1D_x-k+1:dof_1D_x;
IndexSta = 1:k;
Iu = repmat(IndexEnd,k,1);
Iv = repmat(IndexSta,k,1);
GradXFluxJ = GradXFluxJ...
    +sparse([Iv',Iu'],[Iu,Iv],-[Asub,Asup],dof_1D_x,dof_1D_x);


% PGradX denote the flux is taken as \hat{u}=u+
% NGradX denote the flux is taken as \hat{u}=u-
% GradX denote the flux is taken as \hat{u}={u}
PGradX = GradX+GradXFluxC+GradXFluxJ/2;
NGradX = GradX+GradXFluxC-GradXFluxJ/2;
 GradX = GradX+GradXFluxC;

%======================================
% Matrices related to v variable
% vMassV and GradV
%======================================
GradV = (GradX)/2*Lmax/Vmax;

% (vf(v),w(v))_Kv
for Lv=0:nv-1
    xi_v=(( (quad_x+1)/2+Lv)*hv-Vmax); % mapping from [-1,1] to physical domain
    %---------------------------------------------
    % Matrix
    %---------------------------------------------
    % value of local matrix
    val_loc=p_val'*(p_val.*xi_v.*quad_w)*Jacobi_v/2/hv;
    Iu=meshgrid(k*Lv+1:k*(Lv+1));
    vMassV=vMassV+sparse(Iu',Iu,val_loc,dof_1D_v,dof_1D_v);
    
end

end
