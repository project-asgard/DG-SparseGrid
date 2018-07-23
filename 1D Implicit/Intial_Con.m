function [b,EE,BB]=Intial_Con(Lev,Deg,Lmax,pde,FMWT_COMP)
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


% compute (u',v)+1/2*u^{-}[v^{+}-v^{-}]
dof_1D = Deg*nx;

b1=sparse(dof_1D,1);
b2=sparse(dof_1D,1);
b3=sparse(dof_1D,1);

E1=sparse(dof_1D,1);
E2=sparse(dof_1D,1);

B=sparse(dof_1D,1);

%******************************************************
% generate 1D matrix for DG (du/dx,v) by weak form
% Here we only consider the central flux in the scheme
% -(u,dv/dx)+<{u},[v]>
%******************************************************
for LL=0:nx-1

    % current cell
    xi_x=hx*(quad_x/2+1/2+LL); % mapping from [-1,1] to physical domain
    
    f1=pde.f1( xi_x );
    f2=pde.f2( xi_x );
    f3=pde.f3( xi_x );% Added by Lin
    
    
    Iu=[Deg*LL+1:Deg*(LL+1)];
    
    b1(Iu)=p_val'*(quad_w.*f1)*Jacobi_x*sqrt(1/hx)/2;
    b2(Iu)=p_val'*(quad_w.*f2)*Jacobi_x*sqrt(1/hx)/2;
    b3(Iu)=p_val'*(quad_w.*f3)*Jacobi_x*sqrt(1/hx)/2;
    
    
    valE1=pde.e1( xi_x );
    valE2=pde.e2( xi_x );
    
    valB=pde.b( xi_x );

    
    E1(Iu)=p_val'*(quad_w.*valE1)*Jacobi_x*sqrt(1/hx)/2;%*2^(-n)/2*2^(n/2);
    E2(Iu)=p_val'*(quad_w.*valE2)*Jacobi_x*sqrt(1/hx)/2;%*2^(-n)/2*2^(n/2);
   
    B(Iu)=p_val'*(quad_w.*valB)*Jacobi_x*sqrt(1/hx)/2;%*2^(-n)/2*2^(n/2);
    
end

b1=FMWT_COMP*b1;
b2=FMWT_COMP*b2;
b3=FMWT_COMP*b3;


E1=FMWT_COMP*E1;
E2=FMWT_COMP*E2;


B=FMWT_COMP*B;

% structure of b
b = struct('b1',b1,'b2',b2,'b3',b3);

% structure of E
EE = struct('E1',E1,'E2',E2);

% structure of B
BB = struct('B',B);


end