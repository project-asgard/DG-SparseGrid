function EMassX=matrix_coeff_TD(Lev_x,k,Lmin,Lmax,EE,FMWT_blocks)
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

% -----------------------------------------------------------
% note vector quad_x(:) (sampling points) is quad_num by 1,
%      vector quad_w(:) (quadrature weights) is quad_num by 1
% -----------------------------------------------------------
[quad_x,quad_w]=lgwt(quad_num,-1,1);

% --------------------------------
% note p_val(:,:) is quad_num by k
% --------------------------------
p_val = legendre(quad_x,k);

%---------------------------
% Jacobi of variable x and v
% Define Matrices
%---------------------------
nx=2^(Lev_x);
hx=(Lmax-Lmin)/nx;
Jacobi_x=hx;
dof_1D_x=k*nx;

nmax = 2*1024;
use_dense = (dof_1D_x <= nmax);

if (use_dense),
    EMassX=zeros(dof_1D_x,dof_1D_x);
else
    EMassX=sparse(dof_1D_x,dof_1D_x);
end;

%===================================
% Convert E to the DG scaling basis
% So EE has local support
%===================================
%EEold = EE;


if (use_dense),
    EE = full(FMWT_COMP_x)'*full(EE);
else
    EE = FMWT_COMP_x'*EE;
end;

% figure;plot(EE)
for LL=0:nx-1
    i1 = k*LL+1;
    i2 = k*(LL+1);
    m = i2-i1+1;
    
    % --------------------------
    % note Iv(:,:) is   m by m, where m = i2-i1+1
    % Iv = [i1,i1+1, ...,i2;
    %       i1,i1+2, ...,i2;
    %       ...
    %       i1,i1+1, ...,i2]
    % --------------------------
    
    ff=p_val*EE(i1:i2)*sqrt(1/hx);
    
    val=(1/hx)*[p_val'*(quad_w.*p_val.*ff)]*(Jacobi_x/2);
    
    if (use_dense),
        EMassX( i1:i2, i1:i2) = EMassX(i1:i2,i1:i2) + val(1:m,1:m);
    else
        Iv = meshgrid(i1:i2);
        Iu = transpose( Iv );
        EMassX=EMassX+sparse(Iu,Iv,val,dof_1D_x,dof_1D_x);
    end;
end


%***************************************
% Following is for Multiwavelet DG Basis
% Simple way to convert??
%***************************************

% Transfer to multi-DG bases
%
% ------------------------------------------------
% Note   operation is equivalent to
% EMassX = kron( FMWT_COMP_x, FMWT_COMP_x) * EMassX
% ------------------------------------------------

%EMassX=FMWT_COMP_x*EMassX*FMWT_COMP_x';
EMassX = apply_FMWT_blocks(Lev_x, FMWT_blocks, EMassX, left_notrans);
EMassX = apply_FMWT_blocks(Lev_x, FMWT_blocks, EMassX, right_trans);
disp('');



