function EE=PoissonSolve(Lev_x,k,Lmax,rho,DeltaX)
%===============================================================
% Compute EE from Poisson solver
% Input : Lev_x,k,Lmax,rho_0,DeltaX
% Output: EE
%===============================================================

%--Quadrature
quad_num=10;
%---------------

[quad_x,quad_w]=lgwt(quad_num,-1,1);
p_val = legendre(quad_x,k);

%---------------------------
% Jacobi of variable x and v
%---------------------------
nx=2^(Lev_x);hx=Lmax/nx;
Jacobi_x=hx;
dof_1D_x=k*nx;

b_poisson=sparse(2*dof_1D_x,1);
% ff=@(x) 4*pi^2*sin(2*pi*x);
for Lx=0:nx-1
    Iu=[k*Lx+1:k*(Lx+1)];

    ff=1-p_val(:,1:k)*rho(Iu);
    
    b_poisson(dof_1D_x+Iu)=b_poisson(dof_1D_x+Iu)+...
        p_val'*(quad_w.*ff)*Jacobi_x*sqrt(1/hx)/2;

%     xi_x=hx*(quad_x/2+1/2+Lx);
%     b_poisson(dof_1D_x+Iu)=p_val'*(quad_w.*ff(xi_x))*Jacobi_x*sqrt(1/hx)/2;
    
    Meval_x(2*Lx+1:2*Lx+2,Iu)=sqrt(1/hx)*legendre([-1,1],k);
    x_node(2*Lx+1:2*Lx+2)=[Lx*hx,(Lx+1)*hx];
end
% b_poisson(1)=0;

%-----------------
% Handling B.C. 
%-----------------
b_poisson(dof_1D_x+1)=0;
b_poisson(end)=0;

x_poisson=DeltaX\b_poisson;

EE=x_poisson(1:dof_1D_x);

% plot for validating
% figure(11)
% plot(x_node,Meval_x*x_poisson(dof_1D_x+1:end),'r-o');hold on;
% plot(x_node,Meval_x*EE,'b-^');hold off;
% % val=Meval_x*x_poisson(dof_1D_x+1:end);
% % [val(1) val(end)]

%=====================================
% Energy
% (int_x |E|^2 dx)/2
%=====================================
E_energy=0;
for Lx=0:nx-1
    
    fx=sqrt(1/hx)*p_val*EE(k*Lx+1:k*(Lx+1));
    tmp_energy=sum(quad_w.*fx.^2*Jacobi_x/2);
%     E_energy=E_energy+tmp_energy;
end

% E_energy=E_energy/2+Eng.eng;

% % Conservation 
% [Eng.num Eng.mot Eng.entp Eng.eng E_energy]
