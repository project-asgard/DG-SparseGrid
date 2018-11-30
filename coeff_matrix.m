function coeff_matrix(term,domain)

%function [vMassV,GradV,GradX,DeltaX,FluxX,FluxV]=matrix_coeff_TI(Lev_x,Lev_v,k,Lmin,Lmax,Vmin,Vmax,FMWT_COMP_x,FMWT_COMP_v)
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

% Things to add ...

% 1. Choice of flux (may require input C)
% 2. Other BCs
% 3. Picking which term type


%--Quadrature
quad_num=10;
%---------------

% ------------------------------------
% compute the trace values
% note vector p_1[] is 1 by k,
%      vector p_2[] is 1 by k
% ------------------------------------
p_1 = legendre(-1,k);
p_2 = legendre(1,k);

% ----------------------------------
% note quad_x(:) is quad_num by 1
%      quad_w(:) is quad_num by 1
%      p_val(:,:) is quad_num by k
%      Dp_val(:,:) is quad_num by k
% ----------------------------------
[quad_x,quad_w]=lgwt(quad_num,-1,1);
p_val = legendre(quad_x,k);
Dp_val = dlegendre(quad_x,k);

%---------------------------
% Jacobi of variable x and v
% Define Matrices
%---------------------------

nv=2^(Lev_v);
hv=(Vmax-Vmin)/nv;
Jacobi_v=hv;
dof_1D_v=k*nv;

if (use_dense),
    vMassV=zeros(dof_1D_v,dof_1D_v);
    GradV=zeros(dof_1D_v,dof_1D_v);
    FluxV=zeros(dof_1D_v,dof_1D_v);%%
else
    vMassV=sparse(dof_1D_v,dof_1D_v);
    GradV=sparse(dof_1D_v,dof_1D_v);
    FluxV=sparse(dof_1D_v,dof_1D_v);%%
end;



%======================================
% Matrices related to v variable
% vMassV, GradV
%======================================
% (vf(v),w(v))_Kv
vMax=0;
for Lv=0:nv-1
    
    % Get index ranges for ...
    
    % Current element
    % -----------------
    % c=k*Lv+1:k*(Lv+1);
    % -----------------
    c1 = k*Lv+1;
    c2 = k*(Lv+1);
    c = c1:c2;
    
    % Previous element
    % ------------------
    % p=k*(Lv-1)+1:k*Lv;
    % ------------------
    p1 = k*(Lv-1)+1;
    p2 = k*Lv;
    p = p1:p2;
    
    % Later element
    % ----------------------
    % l=k*(Lv+1)+1:k*(Lv+2);
    % ----------------------
    l1 = k*(Lv+1)+1;
    l2 = k*(Lv+2);
    l = l1:l2;
    
    % Map from [-1,1] to physical domain
    % ----------------------------------
    
    xi_v=(( (quad_x+1)/2+Lv)*hv+Vmin);
    
    % Build Average (AVG) and Jump (JMP) operators
    % --------------------------------------------
    
    val_AVG=(1/hv)*[-p_1'*p_2/2  -p_1'*p_1/2,...   % for x1
        p_2'*p_2/2   p_2'*p_1/2];     % for x2
    
    val_JMP=1/hv*[p_1'*p_2 -p_1'*p_1,... % for x1
        -p_2'*p_2  p_2'*p_1]/2; % for x2 %%
    
    
    % 1. Perform volume integral
    % -----------------------
    val_vMassV=p_val'*(p_val.*xi_v.*quad_w)*Jacobi_v/2/hv;
    val_AVG=1/hv*[Dp_val'*(quad_w.*p_val)];
    
    Iu=meshgrid(k*Lv+1:k*(Lv+1));
    
    vMassV=vMassV+sparse(Iu',Iu,val_vMassV,dof_1D_v,dof_1D_v);
    GradV =GradV +sparse(Iu',Iu,val_AVG,dof_1D_v,dof_1D_v);
    
    
    % 2. Setup numerical flux choice (interior elements only)
    % ---------------------------
    
    if Lv<nv-1 && Lv>0
        
        Iu=[meshgrid(p),meshgrid(c),meshgrid(c),meshgrid(l)];
        Iv=[meshgrid(c)',meshgrid(c)',meshgrid(c)',meshgrid(c)'];
        
    end
    
    % 3. Setup boundary conditions
    % -------------------------
    
    BC = 'periodic';
    
    % Periodic BCs
    % ------------
    
    if strcmp(BC,'periodic')
        
        % Element left side
        if Lv==0
            
            Iu=[meshgrid([k*(nv-1)+1:k*(nv)]),meshgrid(c),meshgrid(c),meshgrid(l)];
            Iv=[meshgrid(c)',meshgrid(c)',meshgrid(c)',meshgrid(c)'];
            
            % Element right side
        elseif Lv==nv-1
            
            Iu=[meshgrid(p),meshgrid(c),meshgrid(c),meshgrid([1:k])];
            Iv=[meshgrid(c)',meshgrid(c)',meshgrid(c)',meshgrid(c)'];
        end
        
    end
    
    % Dirichelt BCs
    % -------------
    
    if strcmp(BC,'dirichlet')
    end
    
    % Neumann BCs
    % -----------
    
    if strcmp(BC,'neumann')
    end
    
    % Apply flux choice / BCs
    % -----------------------
    
    GradV=GradV-sparse(Iv,Iu,val_AVG,dof_1D_v,dof_1D_v);
    FluxV=FluxV+sparse(Iv,Iu,val_JMP,dof_1D_v,dof_1D_v);
    
end





%***************************************
% Following is for Multiwavelet DG Basis
% Simple way to convert??
%***************************************
% Transfer to multi-DG bases
vMassV = FMWT_COMP_v*vMassV*FMWT_COMP_v';
GradX = FMWT_COMP_x*GradX*FMWT_COMP_x';
GradV = FMWT_COMP_v*GradV*FMWT_COMP_v';

FluxV = FMWT_COMP_v*FluxV*FMWT_COMP_v';
FluxX = FMWT_COMP_x*FluxX*FMWT_COMP_x';

use_blkdiag = 0;
if (use_blkdiag),
    DeltaX = blkdiag(FMWT_COMP_x,FMWT_COMP_x)*...
        DeltaX*...
        blkdiag(FMWT_COMP_x',FMWT_COMP_x');
else
    % ----------------------------------------
    % note blkdiag(A,B) creates a block diagonal matrix
    % [A, 0;
    % [0, B]
    %
    % let F = FMWT_COMP_x
    %
    % DeltaX = [F, 0;   [D11,  D12;   [F', 0;
    %           0, F] * [D21,  D22] *  0 , F']
    % DeltaX = [ F*D11*F',    F*D12*F';
    %            F*D21*F',    F*D22*F']
    % computed as
    % step 1:  DeltaX = blkdiag(F,F) * DeltaX
    % to form  [F * D11, F * D12;
    %           F * D21, F * D22 ]
    %
    % step 2:  DeltaX = DeltaX * blkdiag(F',F')
    % to form [ (F*D11)*F',   (F*D12)*F';
    %           (F*D21)*F',   (F*D22)*F']
    % ----------------------------------------
    
    % ---------------------------------
    % note shape of DeltaX matrix
    % is  (2*dof_1D_x)  by (2*dof_1D_x)
    % ---------------------------------
    m = dof_1D_x;
    F = FMWT_COMP_x;
    
    % --------------------------------
    % step 1: multiply by blkdiag(F,F) on left
    % --------------------------------
    
    DeltaX(1:m, 1:(2*m)) = F * DeltaX(1:m,1:(2*m));
    DeltaX((m+1):(2*m), 1:(2*m)) = F * DeltaX( (m+1):(2*m), 1:(2*m) );
    
    % -----------------------------------
    % multiply by blkdiag(F',F') on right
    % -----------------------------------
    
    DeltaX(1:(2*m), 1:m) = DeltaX(1:(2*m),1:m) * F';
    DeltaX(1:(2*m), (m+1):(2*m)) = DeltaX(1:(2*m),(m+1):(2*m)) * F';
end;

% vMax

end