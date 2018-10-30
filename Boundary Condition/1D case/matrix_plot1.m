function [Meval_x,x_node]=matrix_plot1(Lev_x,k,Lmax,FMWT_COMP_x)
%=========================================
% Generate the evaluation matrix and plotting points
%=========================================
%-----------------
%--Quadrature
%-----------------
[quad_x,quad_w]=lgwt(k,-1,1);
p_val = legendre(quad_x,k);

%---------------------------
% Jacobi of variable x and v
%---------------------------
nx=2^(Lev_x);hx=Lmax/nx;
dof_1D_x=k*nx;
x_node=zeros(dof_1D_x,1);
Meval_x=sparse(dof_1D_x,dof_1D_x);



for Lx=0:nx-1
    %---------------------------------------------
    % Generate the coefficients for DG bases
    %---------------------------------------------
    Iu=[k*Lx+1:k*(Lx+1)];
    xi_x=hx*(quad_x/2+1/2+Lx);
    
    Meval_x(Iu,Iu)=sqrt(1/hx)*p_val;
    x_node(k*Lx+1:k*Lx+k)=xi_x;

end



%***************************************
% Following is for Multiwavelet DG Basis
% Simple way to convert??
%***************************************




