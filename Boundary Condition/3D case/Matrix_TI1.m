function [GradX,GradY] = Matrix_TI1(Lev,Deg,Lmax,FMWT_COMP_x)
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

[quad_x,quad_w]=lgwt(quad_num,-1,1);
p_val = legendre(quad_x,Deg);
Dp_val = dlegendre(quad_x,Deg);
%---------------------------
% Define Matrices
%---------------------------
nx=2^(Lev);hx=Lmax/nx;
dof_1D_x=Deg*nx;
GradX=sparse(dof_1D_x,dof_1D_x);
GradY=sparse(dof_1D_x,dof_1D_x);
%******************************************************
% generate 1D matrix for DG (du/dx,v) by weak form
% Here we only consider the central flux in the scheme
% -(u,dv/dx)+<{u},[v]>
%******************************************************
%======================================
% Matrices related to x variable GradX
%======================================
% compute (u',v)+1/2*u^{-}[v^{+}-v^{-}]
% generate 1D matrix for DG
for Lx=0:nx-1

    %---------------------------------------------
    % Matrix GradX and EMassX
    %---------------------------------------------
    val=-1/hx*[Dp_val'*(quad_w*ones(1,Deg).*p_val)];
    
    Iu=[meshgrid(Deg*Lx+1:Deg*(Lx+1))]';
    Iv=[meshgrid(Deg*Lx+1:Deg*(Lx+1))];
    GradX=GradX+sparse(Iu,Iv,val,dof_1D_x,dof_1D_x);
    GradY=GradY+sparse(Iu,Iv,val,dof_1D_x,dof_1D_x);
    
    c=Deg*Lx+1:Deg*(Lx+1);
    p=Deg*(Lx-1)+1:Deg*Lx;
    l=Deg*(Lx+1)+1:Deg*(Lx+2);
    p_0=zeros(1,Deg);
    
    val=1/hx*[-p_1'*p_2/2  -p_1'*p_1/2,...   % for x1
               p_2'*p_2/2   p_2'*p_1/2];     % for x2
    
    val1=1/hx*[-p_1'*p_0/2  -p_1'*p_1/2,...   % for x1
                p_2'*p_2/2   p_2'*p_1/2];     % for x2
    
    val2=1/hx*[-p_1'*p_2/2  -p_1'*p_1/2,...   % for x1
                p_2'*p_2/2   p_2'*p_0/2];     % for x2
    
      %computing righthand side & E = cos(pi*x) E(0)=1,E(1)=-1;
    Iv=[meshgrid(c)',meshgrid(c)',meshgrid(c)',meshgrid(c)'];
    
    if Lx<nx-1 && Lx>0
        
        Iu=[meshgrid(p),meshgrid(c),meshgrid(c),meshgrid(l)];
        GradX=GradX+sparse(Iv,Iu,val,dof_1D_x,dof_1D_x);
        
    elseif Lx==0
        
        Iu=[meshgrid(c),meshgrid(c),meshgrid(c),meshgrid(l)];
        GradX=GradX+sparse(Iv,Iu,val1,dof_1D_x,dof_1D_x);
        
    elseif Lx==nx-1
        
        Iu=[meshgrid(p),meshgrid(c),meshgrid(c),meshgrid(c)];
        GradX=GradX+sparse(Iv,Iu,val2,dof_1D_x,dof_1D_x);
    end
    
%     GradX=GradX-sparse(Iv,Iu,val,dof_1D_x,dof_1D_x);
      
    
    
end

GradX=FMWT_COMP_x*GradX*FMWT_COMP_x';
GradY=FMWT_COMP_x*GradY*FMWT_COMP_x';
end