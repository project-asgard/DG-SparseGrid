function [vMassV,GradV,GradX,DeltaX]=matrix_coeff_TI(Lev,k,Lmax,Vmax,FMWT_COMP_x,FMWT_COMP_v)
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

%======================================
% Matrices related to x variable
% GradX and DeltaX
%======================================
% compute (u',v)+1/2*u^{-}[v^{+}-v^{-}]
% generate 1D matrix for DG
for Lx=0:nx-1

    %---------------------------------------------
    % Matrix GradX and EMassX
    %---------------------------------------------
    val=1/hx*[Dp_val'*(quad_w.*p_val)];
    
    Iu=[meshgrid(k*Lx+1:k*(Lx+1))]';
    Iv=[meshgrid(k*Lx+1:k*(Lx+1))];
    GradX=GradX+sparse(Iu,Iv,val,dof_1D_x,dof_1D_x);
    
    DeltaX=DeltaX+sparse([Iu,dof_1D_x+Iu,Iu],[dof_1D_x+Iv,Iv,Iv],...
        [val,val,diag(ones(1,k))],2*dof_1D_x,2*dof_1D_x);
    
    c=k*Lx+1:k*(Lx+1);
    p=k*(Lx-1)+1:k*Lx;
    l=k*(Lx+1)+1:k*(Lx+2);
    
    val=1/hx*[-p_1'*p_2/2  -p_1'*p_1/2,...   % for x1
        p_2'*p_2/2   p_2'*p_1/2];     % for x2
    
    val_u=1/hx*[-p_1'*p_1, p_2'*p_1];
    val_s=1/hx*[-p_1'*p_2, p_2'*p_2];
    
    Iv=[meshgrid(c)',meshgrid(c)',meshgrid(c)',meshgrid(c)'];
    
    if Lx<nx-1 && Lx>0
        
        Iu=[meshgrid(p),meshgrid(c),meshgrid(c),meshgrid(l)];
%         val_u=1/hx*[-p_1'*p_1, p_2'*p_1];
    
    elseif Lx==0
        
        Iu=[meshgrid([k*(nx-1)+1:k*(nx)]),meshgrid(c),meshgrid(c),meshgrid(l)];
%         val_u=1/hx*[-p_1'*p_1, p_2'*p_1];
   
    elseif Lx==nx-1
        
        Iu=[meshgrid(p),meshgrid(c),meshgrid(c),meshgrid([1:k])];
%         val_u=1/hx*[-p_1'*p_1, p_2'*(p_1-p_1)];
    
    end
    
    GradX=GradX-sparse(Iv,Iu,val,dof_1D_x,dof_1D_x);
    DeltaX=DeltaX+sparse([Iv(:,1:2*k),Iv(:,1:2*k)+dof_1D_x],...
        ...[Iu(:,1:2*k)+dof_1D_x,Iu(:,2*k+1:end)],...
        [Iu(:,2*k+1:end)+dof_1D_x,Iu(:,1:2*k)],...
        -[val_u,val_s],2*dof_1D_x,2*dof_1D_x);
    
    
end

% Handle B.C. for Poisson solver
DeltaX(dof_1D_x+1,:)=0;
DeltaX(dof_1D_x+1,dof_1D_x+[1:k])=sqrt(1/hx)*legendre(-1,k);

DeltaX(end,:)=0;
DeltaX(end,end-k+[1:k])=sqrt(1/hx)*legendre(1,k);

%======================================
% Matrices related to v variable
% vMassV and GradV
%======================================
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
    
    val=1/hv*[Dp_val'*(quad_w.*p_val)];
    Iu=[meshgrid(k*Lv+1:k*(Lv+1))]';
    Iv=[meshgrid(k*Lv+1:k*(Lv+1))];
    GradV=GradV+sparse(Iu,Iv,val,dof_1D_v,dof_1D_v);
    
    c=k*Lv+1:k*(Lv+1);
    p=k*(Lv-1)+1:k*Lv;
    l=k*(Lv+1)+1:k*(Lv+2);
    
    val=1/hv*[-p_1'*p_2/2  -p_1'*p_1/2,...   % for x1
        p_2'*p_2/2   p_2'*p_1/2];     % for x2
    
    if Lv<nv-1 && Lv>0
        
        Iu=[meshgrid(p),meshgrid(c),meshgrid(c),meshgrid(l)];
        Iv=[meshgrid(c)',meshgrid(c)',meshgrid(c)',meshgrid(c)'];
        
    elseif Lv==0
        
        Iu=[meshgrid([k*(nv-1)+1:k*(nv)]),meshgrid(c),meshgrid(c),meshgrid(l)];
        Iv=[meshgrid(c)',meshgrid(c)',meshgrid(c)',meshgrid(c)'];
        
    elseif Lv==nv-1
        
        Iu=[meshgrid(p),meshgrid(c),meshgrid(c),meshgrid([1:k])];
        Iv=[meshgrid(c)',meshgrid(c)',meshgrid(c)',meshgrid(c)'];
    end
    GradV=GradV-sparse(Iv,Iu,val,dof_1D_v,dof_1D_v);
end

%***************************************
% Following is for Multiwavelet DG Basis
% Simple way to convert??
%***************************************
% Transfer to multi-DG bases
vMassV = FMWT_COMP_v*vMassV*FMWT_COMP_v';
GradX = FMWT_COMP_x*GradX*FMWT_COMP_x';
GradV = FMWT_COMP_v*GradV*FMWT_COMP_v';

DeltaX = blkdiag(FMWT_COMP_x,FMWT_COMP_x)*...
                DeltaX*...
         blkdiag(FMWT_COMP_x',FMWT_COMP_x');


