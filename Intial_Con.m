function [f_v,f_x]=Intial_Con(Lev_x,Lev_v,k,Lmax,Vmax,pde,...
    FMWT_COMP_x,FMWT_COMP_v)
%==================================================================
% This code computes the initial conditions for f
% Input:
%   Lev_x,Lev_v,k,Lmax,Vmax,pde 
% Output: 
%   f_v, f_x
%==================================================================


% Get the Legendre-Gauss nodes (quad_x) and weights (quad_w) on the domain
% [-1,+1] for performing quadrature.

quad_num = 10;
[quad_x,quad_w] = lgwt(quad_num,-1,1);

% Get the Legendre basis function evaluated at the Legendre-Gauss nodes up
% to order k.

p_val = transpose(legendre(quad_x,k));

% Get grid spacing for both  x and v.

nx=2^(Lev_x);
hx=Lmax/nx;
dof_1D_x=k*nx;
f_x=zeros(dof_1D_x,1);

nv=2^(Lev_v);
hv=2*Vmax/nv;
dof_1D_v=k*nv;
f_v=zeros(dof_1D_v,1);

% Initial Condition for f_x
params = pde.params;

for i=0:nx-1
    
    % Map quad_x from [-1,+1] to [0,LMax] physical domain.
    xi_x = hx*(quad_x/2+1/2+i);
    
    % Get the f(x) initial condition at the quadrature points.
    fxHere = pde.Fx_0(xi_x, params);
       
    % Generate the coefficients for DG basis    
    for thisk=1:k
        
        this_k_legendre = p_val(thisk,:);
        this_quad = (quad_w .* fxHere);
        f_x(k*i+thisk) = mtimes(this_k_legendre,this_quad) * hx * sqrt(1/hx)/2;
        
    end
    
end

% Initial Condition for f_v

for i=0:nv-1
    
    % Map fron quad_v from [-1,1] to [-Vmax,+Vmax] domain.
    xi_v=hv*(quad_x/2+1/2+i)-Vmax;
    
    % Get the f(v) initial condition at the quadrature points
    fvHere = pde.Fv_0(xi_v, params);
    
    % Generate the coefficients for DG basis 
    for thisk=1:k
        
        this_k_legendre = p_val(thisk,:);
        this_quad = (quad_w .* fvHere);
        f_v(k*i+thisk) = mtimes(this_k_legendre,this_quad) * hv * sqrt(1/hv)/2;
        
    end   
       
end


% Transfer to multi-DG bases

f_v = mtimes( FMWT_COMP_v, f_v );
f_x = mtimes( FMWT_COMP_x, f_x );


