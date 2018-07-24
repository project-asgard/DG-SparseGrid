function [f_v,f_x]=Intial_Con(Lev_x,Lev_v,k,Lmax,Vmax,pde,...
    FMWT_COMP_x,FMWT_COMP_v)
%==================================================================
% This code computes the initial conditions for f and rho=int_v fdv
% Input: Lev_x,Lev_v,k,Lmax,Vmax,pde
% Output: f_v, f_x, rho, and Eng
%
%==================================================================
%--Quadrature
quad_num=10;
%---------------
[quad_x,quad_w]=lgwt(quad_num,-1,1);

p_val = legendre(quad_x,k);

%---------------------------
% Jacobi of variable x and v
%---------------------------
nx=2^(Lev_x);hx=Lmax/nx;
Jacobi_x=hx;
dof_1D_x=k*nx;
f_x=zeros(dof_1D_x,1);

nv=2^(Lev_v);hv=2*Vmax/nv;
Jacobi_v=hv;
dof_1D_v=k*nv;
f_v=zeros(dof_1D_v,1);

%=================================================
% Initial Condition for f_x
%=================================================
for Lx=0:nx-1
    xi_x=hx*(quad_x/2+1/2+Lx); % mapping from [-1,1] to physical domain
    %---------------------------------------------
    % Generate the coefficients for DG bases
    %---------------------------------------------
    Iu=[k*Lx+1:k*(Lx+1)];
    
    f_x(Iu)=p_val'*(quad_w.*pde.Fx_0(xi_x))*Jacobi_x*sqrt(1/hx)/2;
end

%=================================================
% Initial Condition for f_v
%=================================================
for Lv=0:nv-1
    xi_v=(( (quad_x+1)/2+Lv)*hv-Vmax); % mapping from [-1,1] to physical domain
    %---------------------------------------------
    % Generate the coefficients for DG bases
    %---------------------------------------------
    Iv=[k*Lv+1:k*(Lv+1)];
    
    f_v(Iv)=p_val'*(quad_w.*pde.Fv_0(xi_v))*Jacobi_v*sqrt(1/hv)/2;
end


%***************************************
% Following is for Multiwavelet DG Basis
%***************************************

% Transfer to multi-DG bases
f_v=FMWT_COMP_v*f_v;
f_x=FMWT_COMP_x*f_x;


