function EMassX=matrix_coeff_TD(Lev_x,k,Lmax,EE,FMWT_COMP_x)
%=============================================================
% Generate time-dependent coefficient matrix
% Vlasolv Solver:   
%   Operators:  EMassX: int_x E(x)*m_i(x)*m_j(x)dx
% Input: Lev, k, dim, Lmax, Vmax
% Note: E is solved by Poisson or Maxwell's equation
%=============================================================

%--Quadrature
quad_num=10;
%---------------

[quad_x,quad_w]=lgwt(quad_num,-1,1);
p_val = legendre(quad_x,k);

%---------------------------
% Jacobi of variable x and v
% Define Matrices
%---------------------------
nx=2^(Lev_x);hx=Lmax/nx;
Jacobi_x=hx;
dof_1D_x=k*nx;

EMassX=sparse(dof_1D_x,dof_1D_x);

%===================================
% Convert E to the DG scaling basis
% So EE has local support
%===================================
EE = FMWT_COMP_x'*EE;

% figure;plot(EE)
for LL=0:nx-1
    Iu=[meshgrid(k*LL+1:k*(LL+1))]';
    Iv=[meshgrid(k*LL+1:k*(LL+1))];

    % Generate EE
%     ff=legendre(0,k)*EE(k*LL+1:k*(LL+1)); % use the middle point-value
    ff=p_val*EE(k*LL+1:k*(LL+1));
    
    val=1/hx*[p_val'*(quad_w.*p_val.*ff)]*Jacobi_x/2;
    
    EMassX=EMassX+sparse(Iu,Iv,val,dof_1D_x,dof_1D_x);
end

%***************************************
% Following is for Multiwavelet DG Basis
% Simple way to convert??
%***************************************
    
% Transfer to multi-DG bases
EMassX=FMWT_COMP_x*EMassX*FMWT_COMP_x';
   


