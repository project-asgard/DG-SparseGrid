function [b,EE]=Intial_Con1(Lev,Deg,Lmax,pde,FMWT_COMP_x)
%========================================================
% Generate the Initial Condition for Maxwells' Eq 
% Input: 
%   Lev denotes the Level for the mesh
%   Deg denotes the degree for polynomial
%   pde gives the rhs and boundary conditions
% Output:
%   L_MWDG: Grad Matrix
% PDE dependent Stiffness Matrix
%========================================================

%--DG parameters
quad_num=10;
%---------------

[quad_x,quad_w]=lgwt(quad_num,-1,1);
p_val = legendre(quad_x,Deg);

nx=2^(Lev);hx=Lmax/nx;
Jacobi_x=hx;

p_1 = legendre(-1,Deg);
p_2 = legendre(1,Deg);
s_1 = pde.E(0)*p_1/2*sqrt(1/hx);
s_2 = -pde.E(1)*p_2/2*sqrt(1/hx);

% compute (u',v)+1/2*u^{-}[v^{+}-v^{-}]
dof_1D = Deg*nx;

b=sparse(dof_1D,1);

E=sparse(dof_1D,1);

%******************************************************
% generate 1D matrix for DG (du/dx,v) by weak form
% Here we only consider the central flux in the scheme
% -(u,dv/dx)+<{u},[v]>
%******************************************************
for LL=0:nx-1

    % current cell
    xi_x=hx*(quad_x/2+1/2+LL); % mapping from [-1,1] to physical domain
    
    f=pde.f( xi_x );

    Iu=[Deg*LL+1:Deg*(LL+1)];
    
    b(Iu)=p_val'*(quad_w.*f)*Jacobi_x*sqrt(1/hx)/2;
   
    
    valE=pde.E( xi_x );

    E(Iu)=p_val'*(quad_w.*valE)*Jacobi_x*sqrt(1/hx)/2;%*2^(-n)/2*2^(n/2);
   

end

b(1:Deg,1)=b(1:Deg,1)+s_1';
b(Deg*(nx-1)+1:Deg*nx,1)=b(Deg*(nx-1)+1:Deg*nx,1)+s_2';


b=FMWT_COMP_x*b;
E=FMWT_COMP_x*E;
% structure of b
b = struct('b',b);

% structure of E
EE = struct('E',E);




end