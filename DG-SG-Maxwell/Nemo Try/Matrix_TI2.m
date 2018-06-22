function GradX = Matrix_TI2(Lev,Deg,Lmax,FMWT_COMP_x)
%========================================================
% Construct the matrix for curl operator on [0,Lmax]
% Input: 
%   Lev denotes the Level for the mesh
%   Deg denotes the degree for polynomial
%   Lmax denotes the domain for x variable
%   FMWT_COMP_x denotes the converting relationship
% Output:
%   GradX: Grad Matrix
%========================================================
%--DG parameters
quad_num=10;
%---------------


% compute the trace values
p_1 = legendre(-1,Deg);
p_2 = legendre(1,Deg);
p_0 = legendre(0,Deg);

[quad_x,quad_w]=lgwt(quad_num,-1,1);
p_val = legendre(quad_x,Deg);
Dp_val = dlegendre(quad_x,Deg);

%---------------------------
% Define Matrices
%---------------------------
nx=2^(Lev);hx=Lmax/nx;
dof_1D_x=Deg*nx;
GradX=sparse(dof_1D_x,dof_1D_x);

%******************************************************
% generate 1D matrix for DG (du/dx,v) by weak form
% Here we only consider the Lax-Friedrichs flux in the scheme
% -(u,dv/dx)+<{u},[v]>
%******************************************************
%======================================
% Matrices related to x variable GradX
%======================================
% generate 1D matrix for DG
for Lx=0:nx-1 

    %---------------------------------------------
    % Matrix GradX and EMassX
    %---------------------------------------------
    val=-1/hx*[Dp_val'*(quad_w.*p_val)];
    
    Iu=[meshgrid(Deg*Lx+1:Deg*(Lx+1))]';
    Iv=[meshgrid(Deg*Lx+1:Deg*(Lx+1))];
    GradX=GradX+sparse(Iu,Iv,val,dof_1D_x,dof_1D_x);
    
    c=Deg*Lx+1:Deg*(Lx+1);
    p=Deg*(Lx-1)+1:Deg*Lx;
    l=Deg*(Lx+1)+1:Deg*(Lx+2);
    
    val=1/hx*[-p_1'*p_2,p_2'*p_2];     

    
    Iv=[meshgrid(c)',meshgrid(c)'];
    
    if Lx<nx-1 && Lx>0
        
        Iu=[meshgrid(p),meshgrid(c)];
    
    elseif Lx==0
        
        Iu=[meshgrid([Deg*(nx-1)+1:Deg*(nx)]),meshgrid(c)];
   
    elseif Lx==nx-1
        
        Iu=[meshgrid(p),meshgrid(c)];
    
    end
    
%     GradX=GradX-sparse(Iv,Iu,val,dof_1D_x,dof_1D_x);
    GradX=GradX+sparse(Iv,Iu,val,dof_1D_x,dof_1D_x);
    
    
end


GradX = FMWT_COMP_x*GradX*FMWT_COMP_x';


end