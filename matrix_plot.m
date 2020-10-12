function [Meval_v,v_node,Meval_x,x_node]=matrix_plot(Lev_x,Lev_v,k,Lmin,Lmax,Vmin,Vmax,FMWT_blocks)
%=========================================
% Generate the evaluation matrix and plotting points
%=========================================
%-----------------
%--Quadrature
%-----------------
[quad_x,quad_w]=lgwt(k,-1,1);
p_val = lin_legendre(quad_x,k);

%---------------------------
% Jacobi of variable x and v
%---------------------------
nx=2^(Lev_x);hx=(Lmax-Lmin)/nx;
dof_1D_x=k*nx;
x_node=zeros(dof_1D_x,1);
Meval_x=sparse(dof_1D_x,dof_1D_x);

nv=2^(Lev_v);hv=(Vmax-Vmin)/nv;
dof_1D_v=k*nv;
v_node=zeros(dof_1D_v,1);
Meval_v=sparse(dof_1D_v,dof_1D_v);

for Lx=0:nx-1
    %---------------------------------------------
    % Generate the coefficients for DG bases
    %---------------------------------------------
    Iu=[k*Lx+1:k*(Lx+1)];
    xi_x=(quad_x/2+1/2+Lx)*hx+Lmin;
    
    Meval_x(Iu,Iu)=sqrt(1/hx)*p_val;
    x_node(k*Lx+1:k*Lx+k)=xi_x;

end

for Lv=0:nv-1
    xi_v=((quad_x+1)/2+Lv)*hv+Vmin; % mapping from [-1,1] to physical domain
    %---------------------------------------------
    % Coefficients for DG bases
    %---------------------------------------------
    Iu=[k*Lv+1:k*(Lv+1)];
    
    Meval_v(Iu,Iu)=sqrt(1/hv)*p_val;
    v_node(Iu)=xi_v;
end

%***************************************
% Following is for Multiwavelet DG Basis
% Simple way to convert??
%***************************************

right_trans = 'RT';
Meval_v = apply_FMWT_blocks(Lev_v, FMWT_blocks, Meval_v, right_trans);
Meval_x = apply_FMWT_blocks(Lev_x, FMWT_blocks, Meval_x, right_trans);


