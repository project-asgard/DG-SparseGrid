function [EE,u]=PoissonSolve(Lev_x,k,Lmax,f,DeltaX,FMWT_COMP_x,Vmax)
%===============================================================
% Compute EE from Poisson solver
% Input : Lev_x,k,Lmax,rho_0,DeltaX
% Output: EE
%===============================================================

%---------------------------
% Jacobi of variable x and v
%---------------------------
nx=2^(Lev_x);hx=Lmax/nx;
dof_1D_x=k*nx;

b_poisson=sparse(2*dof_1D_x,1);


%-----------------------------------------------------
% Handling B.C. for enforcing phi(0)=phi(Lmax)=0
%-----------------------------------------------------
% Phi_v(0,0) corresponding to the

% f is the given coeffients of the wavelet basis
% rho = int_v f dv: 
% we use the orthogonal property (1,w_j) = 0 for the j>1
%       and (1,w_1) = sqrt(2*Vmax)
% RHS = int_x (1-rho)*q dx 
% by orthogonal property, only the first basis gives nonzero RHS, which is
%       sqrt(Lmax)


Index=[1:k]'+[0:k^2:(2^(Lev_x)-1)*k^2]; 
% take out all the coefs corresponding to the first basis corresponding to v
% compute (1-rho,q)
tmp_b=[sqrt(Lmax); zeros(dof_1D_x-1,1)]...
      -sqrt(2*Vmax)*f(Index(:));%f(1:dof_1D_x);

% Boundary condition is enforced as phi(0) = 0; phi(Lmax) = 0
% 1. Convert back to Lengendre basis: C1 = FMWT_COMP_x(:,2:end-1)')*tmp_b
% We exclude first and last basis, and assume the coefs of Lengendre basis are zeros
% 2. Convert to wavelet basis: C2 = FMWT_COMP_x(:,2:end-1)*C1
%b_poisson(dof_1D_x+1:end)=(FMWT_COMP_x(:,2:end-1)*FMWT_COMP_x(:,2:end-1)')*tmp_b;
% Suggested by Ed, the above is changed to
b_poisson(dof_1D_x+1:end)=FMWT_COMP_x(:,2:end-1)*(FMWT_COMP_x(:,2:end-1)'*tmp_b);


x_poisson=DeltaX\b_poisson;

EE=x_poisson(1:dof_1D_x);



u = x_poisson(1+dof_1D_x:end);
