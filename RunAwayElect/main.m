% main code for LDG Poisson equation
format short e
addpath(genpath(pwd))

Lev = 4;
Deg = 2;

Lstart = 0;
Lend = 1;
Lmax = Lend-Lstart;

func = @(x)(sin(x));

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
% GradX=sparse(dof_1D_x,dof_1D_x);
DeltaX=sparse(2*dof_1D_x,2*dof_1D_x);

b_poisson=sparse(2*dof_1D_x,1);

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
    
    Iu=[meshgrid(Deg*Lx+1:Deg*(Lx+1))]';
    Iv=[meshgrid(Deg*Lx+1:Deg*(Lx+1))];
%     GradX=GradX+sparse(Iu,Iv,val,dof_1D_x,dof_1D_x);
    
    DeltaX=DeltaX+sparse([Iu,dof_1D_x+Iu,Iu],[dof_1D_x+Iv,Iv,Iv],...
        [val,val,diag(ones(1,Deg))],2*dof_1D_x,2*dof_1D_x);
    
    c=Deg*Lx+1:Deg*(Lx+1);
    p=Deg*(Lx-1)+1:Deg*Lx;
    l=Deg*(Lx+1)+1:Deg*(Lx+2);
    
    val=1/hx*[-p_1'*p_2/2  -p_1'*p_1/2,...   % for x1
        p_2'*p_2/2   p_2'*p_1/2];     % for x2
    
    val_u=1/hx*[-p_1'*p_1, p_2'*p_1];
    val_s=1/hx*[-p_1'*p_2, p_2'*p_2];
    
    Iv=[meshgrid(c)',meshgrid(c)',meshgrid(c)',meshgrid(c)'];
    
    if Lx<nx-1 && Lx>0
        
        Iu=[meshgrid(p),meshgrid(c),meshgrid(c),meshgrid(l)];
%         val_u=1/hx*[-p_1'*p_1, p_2'*p_1];
    
    elseif Lx==0
        
        Iu=[meshgrid([Deg*(nx-1)+1:Deg*(nx)]),meshgrid(c),meshgrid(c),meshgrid(l)];
%         val_u=1/hx*[-p_1'*p_1, p_2'*p_1];
   
    elseif Lx==nx-1
        
        Iu=[meshgrid(p),meshgrid(c),meshgrid(c),meshgrid([1:Deg])];
%         val_u=1/hx*[-p_1'*p_1, p_2'*(p_1-p_1)];
    
    end
    
%     GradX=GradX-sparse(Iv,Iu,val,dof_1D_x,dof_1D_x);
    DeltaX=DeltaX+sparse([Iv(:,1:2*Deg),Iv(:,1:2*Deg)+dof_1D_x],...
        ...[Iu(:,1:2*k)+dof_1D_x,Iu(:,2*k+1:end)],...
        [Iu(:,2*Deg+1:end)+dof_1D_x,Iu(:,1:2*Deg)],...
        -[val_u,val_s],2*dof_1D_x,2*dof_1D_x);
    
%     % get xhat
%     MidPoint = ((Lstart+(Lx+1)*hx)+(Lstart+(Lx)*hx) )/2;
%     xhat = (x-MidPoint)*2/hx;
%     
%     val = func(xhat);
%     b_poisson = b_poisson+sparse(Iv(:,1:2*Deg)+dof_1D_x,ones(Deg,1),val,2*dof_1D_x,1);
    
end

% Handle B.C. for Poisson solver
DeltaX(dof_1D_x+1,:)=0;
DeltaX(dof_1D_x+1,dof_1D_x+[1:Deg])=sqrt(1/hx)*legendre(-1,Deg);

DeltaX(end,:)=0;
DeltaX(end,end-Deg+[1:Deg])=sqrt(1/hx)*legendre(1,Deg);

A11 = DeltaX(1:dof_1D_x,1:dof_1D_x);
A12 = DeltaX(1:dof_1D_x,dof_1D_x+1:end);
A21 = DeltaX(dof_1D_x+1:end,1:dof_1D_x);
A22 = DeltaX(dof_1D_x+1:end,1+dof_1D_x:end);

A12*A21

% backward Euler



