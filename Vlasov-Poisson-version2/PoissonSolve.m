function EE=PoissonSolve(Lev_x,k,Lmax,f,DeltaX,FMWT_COMP_x,Vmax)
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

Index=[1:k]'+[0:k^2:(2^(Lev_x)-1)*k^2];
tmp_b=[sqrt(Lmax); zeros(dof_1D_x-1,1)]...
      -sqrt(2*Vmax)*f(Index(:));%f(1:dof_1D_x);


b_poisson(dof_1D_x+1:end)=(FMWT_COMP_x(:,2:end-1)*FMWT_COMP_x(:,2:end-1)')*tmp_b;


x_poisson=DeltaX\b_poisson;

EE=x_poisson(1:dof_1D_x);

