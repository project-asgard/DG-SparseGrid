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
b = sparse(dof_1D,1);

b3x1=sparse(dof_1D,1);
b3x2=sparse(dof_1D,1);
b3x3=sparse(dof_1D,1);
b3y1=sparse(dof_1D,1);
b3y2=sparse(dof_1D,1);
b3y3=sparse(dof_1D,1);
b3z1=sparse(dof_1D,1);
b3z2=sparse(dof_1D,1);
b3z3=sparse(dof_1D,1);

b4x1=sparse(dof_1D,1);
b4x2=sparse(dof_1D,1);
b4x3=sparse(dof_1D,1);
b4y1=sparse(dof_1D,1);
b4y2=sparse(dof_1D,1);
b4y3=sparse(dof_1D,1);
b4z1=sparse(dof_1D,1);
b4z2=sparse(dof_1D,1);
b4z3=sparse(dof_1D,1);

Ex1=sparse(dof_1D,1);
Ex2=sparse(dof_1D,1);
Ex3=sparse(dof_1D,1);
Ey1=sparse(dof_1D,1);
Ey2=sparse(dof_1D,1);
Ey3=sparse(dof_1D,1);
Ez1=sparse(dof_1D,1);
Ez2=sparse(dof_1D,1);
Ez3=sparse(dof_1D,1);

Bx1=sparse(dof_1D,1);
Bx2=sparse(dof_1D,1);
Bx3=sparse(dof_1D,1);
By1=sparse(dof_1D,1);
By2=sparse(dof_1D,1);
By3=sparse(dof_1D,1);
Bz1=sparse(dof_1D,1);
Bz2=sparse(dof_1D,1);
Bz3=sparse(dof_1D,1);

%******************************************************
% generate 1D matrix for DG (du/dx,v) by weak form
% Here we only consider the central flux in the scheme
% -(u,dv/dx)+<{u},[v]>
%******************************************************
for LL=0:nx-1

    % current cell
    xi_x=hx*(quad_x/2+1/2+LL); % mapping from [-1,1] to physical domain
    
    f3x=pde.rhs3x( xi_x );
    f3y=pde.rhs3y( xi_x );
    f3z=pde.rhs3z( xi_x );
    
    f4x=pde.rhs4x( xi_x );
    f4y=pde.rhs4y( xi_x );
    f4z=pde.rhs4z( xi_x );
    
    Iu=[Deg*LL+1:Deg*(LL+1)];
    
    b3x1(Iu)=p_val'*(quad_w.*f3x(:,1))*Jacobi_x*sqrt(1/hx)/2;
    b3x2(Iu)=p_val'*(quad_w.*f3x(:,2))*Jacobi_x*sqrt(1/hx)/2;
    b3x3(Iu)=p_val'*(quad_w.*f3x(:,3))*Jacobi_x*sqrt(1/hx)/2;
    b3y1(Iu)=p_val'*(quad_w.*f3y(:,1))*Jacobi_x*sqrt(1/hx)/2;
    b3y2(Iu)=p_val'*(quad_w.*f3y(:,2))*Jacobi_x*sqrt(1/hx)/2;
    b3y3(Iu)=p_val'*(quad_w.*f3y(:,3))*Jacobi_x*sqrt(1/hx)/2;
    b3z1(Iu)=p_val'*(quad_w.*f3z(:,1))*Jacobi_x*sqrt(1/hx)/2;
    b3z2(Iu)=p_val'*(quad_w.*f3z(:,2))*Jacobi_x*sqrt(1/hx)/2;
    b3z3(Iu)=p_val'*(quad_w.*f3z(:,3))*Jacobi_x*sqrt(1/hx)/2;
    
    b4x1(Iu)=p_val'*(quad_w.*f4x(:,1))*Jacobi_x*sqrt(1/hx)/2;
    b4x2(Iu)=p_val'*(quad_w.*f4x(:,2))*Jacobi_x*sqrt(1/hx)/2;
    b4x3(Iu)=p_val'*(quad_w.*f4x(:,3))*Jacobi_x*sqrt(1/hx)/2;
    b4y1(Iu)=p_val'*(quad_w.*f4y(:,1))*Jacobi_x*sqrt(1/hx)/2;
    b4y2(Iu)=p_val'*(quad_w.*f4y(:,2))*Jacobi_x*sqrt(1/hx)/2;
    b4y3(Iu)=p_val'*(quad_w.*f4y(:,3))*Jacobi_x*sqrt(1/hx)/2;
    b4z1(Iu)=p_val'*(quad_w.*f4z(:,1))*Jacobi_x*sqrt(1/hx)/2;
    b4z2(Iu)=p_val'*(quad_w.*f4z(:,2))*Jacobi_x*sqrt(1/hx)/2;
    b4z3(Iu)=p_val'*(quad_w.*f4z(:,3))*Jacobi_x*sqrt(1/hx)/2;
   
    
    valEx=pde.Ex( xi_x );
    valEy=pde.Ey( xi_x );
    valEz=pde.Ez( xi_x );
    
    valBx=pde.Bx( xi_x );
    valBy=pde.By( xi_x );
    valBz=pde.Bz( xi_x );
    
    
    Ex1(Iu)=p_val'*(quad_w.*valEx(:,1))*Jacobi_x*sqrt(1/hx)/2;%*2^(-n)/2*2^(n/2);
    Ex2(Iu)=p_val'*(quad_w.*valEx(:,2))*Jacobi_x*sqrt(1/hx)/2;%*2^(-n)/2*2^(n/2);
    Ex3(Iu)=p_val'*(quad_w.*valEx(:,3))*Jacobi_x*sqrt(1/hx)/2;%*2^(-n)/2*2^(n/2);
    Ey1(Iu)=p_val'*(quad_w.*valEy(:,1))*Jacobi_x*sqrt(1/hx)/2;%*2^(-n)/2*2^(n/2);
    Ey2(Iu)=p_val'*(quad_w.*valEy(:,2))*Jacobi_x*sqrt(1/hx)/2;%*2^(-n)/2*2^(n/2);
    Ey3(Iu)=p_val'*(quad_w.*valEy(:,3))*Jacobi_x*sqrt(1/hx)/2;%*2^(-n)/2*2^(n/2);
    Ez1(Iu)=p_val'*(quad_w.*valEz(:,1))*Jacobi_x*sqrt(1/hx)/2;%*2^(-n)/2*2^(n/2);
    Ez2(Iu)=p_val'*(quad_w.*valEz(:,2))*Jacobi_x*sqrt(1/hx)/2;%*2^(-n)/2*2^(n/2);
    Ez3(Iu)=p_val'*(quad_w.*valEz(:,3))*Jacobi_x*sqrt(1/hx)/2;%*2^(-n)/2*2^(n/2);
    
    Bx1(Iu)=p_val'*(quad_w.*valBx(:,1))*Jacobi_x*sqrt(1/hx)/2;%*2^(-n)/2*2^(n/2);
    Bx2(Iu)=p_val'*(quad_w.*valBx(:,2))*Jacobi_x*sqrt(1/hx)/2;%*2^(-n)/2*2^(n/2);
    Bx3(Iu)=p_val'*(quad_w.*valBx(:,3))*Jacobi_x*sqrt(1/hx)/2;%*2^(-n)/2*2^(n/2);
    By1(Iu)=p_val'*(quad_w.*valBy(:,1))*Jacobi_x*sqrt(1/hx)/2;%*2^(-n)/2*2^(n/2);
    By2(Iu)=p_val'*(quad_w.*valBy(:,2))*Jacobi_x*sqrt(1/hx)/2;%*2^(-n)/2*2^(n/2);
    By3(Iu)=p_val'*(quad_w.*valBy(:,3))*Jacobi_x*sqrt(1/hx)/2;%*2^(-n)/2*2^(n/2);
    Bz1(Iu)=p_val'*(quad_w.*valBz(:,1))*Jacobi_x*sqrt(1/hx)/2;%*2^(-n)/2*2^(n/2);
    Bz2(Iu)=p_val'*(quad_w.*valBz(:,2))*Jacobi_x*sqrt(1/hx)/2;%*2^(-n)/2*2^(n/2);
    Bz3(Iu)=p_val'*(quad_w.*valBz(:,3))*Jacobi_x*sqrt(1/hx)/2;%*2^(-n)/2*2^(n/2);
    

end

b3x1=FMWT_COMP*b3x1;
b3x2=FMWT_COMP*b3x2;
b3x3=FMWT_COMP*b3x3;
b3y1=FMWT_COMP*b3y1;
b3y2=FMWT_COMP*b3y2;
b3y3=FMWT_COMP*b3y3;
b3z1=FMWT_COMP*b3z1;
b3z2=FMWT_COMP*b3z2;
b3z3=FMWT_COMP*b3z3;

b4x1=FMWT_COMP*b4x1;
b4x2=FMWT_COMP*b4x2;
b4x3=FMWT_COMP*b4x3;
b4y1=FMWT_COMP*b4y1;
b4y2=FMWT_COMP*b4y2;
b4y3=FMWT_COMP*b4y3;
b4z1=FMWT_COMP*b4z1;
b4z2=FMWT_COMP*b4z2;
b4z3=FMWT_COMP*b4z3;

Ex1=FMWT_COMP*Ex1;
Ex2=FMWT_COMP*Ex2;
Ex3=FMWT_COMP*Ex3;
Ey1=FMWT_COMP*Ey1;
Ey2=FMWT_COMP*Ey2;
Ey3=FMWT_COMP*Ey3;
Ez1=FMWT_COMP*Ez1;
Ez2=FMWT_COMP*Ez2;
Ez3=FMWT_COMP*Ez3;

Bx1=FMWT_COMP*Bx1;
Bx2=FMWT_COMP*Bx2;
Bx3=FMWT_COMP*Bx3;
By1=FMWT_COMP*By1;
By2=FMWT_COMP*By2;
By3=FMWT_COMP*By3;
Bz1=FMWT_COMP*Bz1;
Bz2=FMWT_COMP*Bz2;
Bz3=FMWT_COMP*Bz3;

% structure of b
b = struct('b3x1',b3x1,'b3x2',b3x2,'b3x3',b3x3,...
    'b3y1',b3y1,'b3y2',b3y2,'b3y3',b3y3,...
    'b3z1',b3z1,'b3z2',b3z2,'b3z3',b3z3,...
    'b4x1',b4x1,'b4x2',b4x2,'b4x3',b4x3,...
    'b4y1',b4y1,'b4y2',b4y2,'b4y3',b4y3,...
    'b4z1',b4z1,'b4z2',b4z2,'b4z3',b4z3);

% structure of E
EE = struct('Ex1',Ex1,'Ex2',Ex2,'Ex3',Ex3,...
    'Ey1',Ey1,'Ey2',Ey2,'Ey3',Ey3,...
    'Ez1',Ez1,'Ez2',Ez2,'Ez3',Ez3);

% structure of B
BB = struct('Bx1',Bx1,'Bx2',Bx2,'Bx3',Bx3,...
    'By1',By1,'By2',By2,'By3',By3,...
    'Bz1',Bz1,'Bz2',Bz2,'Bz3',Bz3);



end