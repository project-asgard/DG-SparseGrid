function [f_v,f_x,rho,Eng]=Intial_Con(Lev_x,Lev_v,k,Lmax,Vmax,pde,...
    FMWT_COMP_x,FMWT_COMP_v)
%==================================================================
% This code computes the initial conditions for f and rho=int_v fdv
% Input: Lev_x,Lev_v,k,Lmax,Vmax,pde
% Output: f_v, f_x, rho, and Eng
%
%==================================================================
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
f_x=zeros(dof_1D_x,1);

nv=2^(Lev_v);hv=2*Vmax/nv;
Jacobi_v=hv;
dof_1D_v=k*nv;
f_v=zeros(dof_1D_v,1);

%=================================================
% Initial Condition for f_x
%=================================================
for Lx=0:nx-1
    xi_x=hx*(quad_x/2+1/2+Lx); % mapping from [-1,1] to physical domain
    %---------------------------------------------
    % Generate the coefficients for DG bases
    %---------------------------------------------
    Iu=[k*Lx+1:k*(Lx+1)];
    
    f_x(Iu)=p_val'*(quad_w.*pde.Fx_0(xi_x))*Jacobi_x*sqrt(1/hx)/2;
end

%=================================================
% Initial Condition for f_v
%=================================================
for Lv=0:nv-1
    xi_v=(( (quad_x+1)/2+Lv)*hv-Vmax); % mapping from [-1,1] to physical domain
    %---------------------------------------------
    % Generate the coefficients for DG bases
    %---------------------------------------------
    Iv=[k*Lv+1:k*(Lv+1)];
    
    f_v(Iv)=p_val'*(quad_w.*pde.Fv_0(xi_v))*Jacobi_v*sqrt(1/hv)/2;
end

%=================================================
% Compute rho and do the following
% check the energy
% 1. partical number int_x rho dx
% 2. momentum int_x int_v vf dvdx
% 3. entropy int_x int_v |f|^2 dvdx
% 4. Total Energy (int_x int_v f|v|^2 dvdx)/2+(int_x |E|^2 dx)/2
%=================================================
rho=zeros(dof_1D_x,1);

E_momentum=0;E_entropy=0;E_energy=0;Num_part=0;
for Lx=0:nx-1
    
    for kx=1:k
        Iu=k*Lx+kx;
        fx=f_x(Iu)*sqrt(1/hx);% for every basis function kx
        
        for Lv=0:nv-1
            xi_v=(( (quad_x+1)/2+Lv)*hv-Vmax);
            for kv=1:k
                Iv=k*Lv+kv;
                fv=sqrt(1/hv)*f_v(Iv)*p_val(:,kv)*fx;
                
                pos=Iu;
                tmp_rho=sum(quad_w.*fv*Jacobi_v/2);
                rho(pos)=rho(pos)+tmp_rho;
                
                %                 xi_v=(( (quad_x+1)/2+Lv)*hv-Vmax);
                
                tmp_momentum=sum(quad_w.*xi_v.*fv*Jacobi_v/2)*sum(quad_w.*p_val(:,kx)*Jacobi_x/2);
                tmp_entropy=sum(quad_w.*fv.^2*Jacobi_v/2)*sum(quad_w.*p_val(:,kx)*Jacobi_x/2);
                tmp_energy=sum(quad_w.*fv.*xi_v.^2*Jacobi_v/2)*sum(quad_w.*p_val(:,kx)*Jacobi_x/2);
                
                E_momentum=E_momentum+tmp_momentum;
                E_entropy=E_entropy+tmp_entropy;
                E_energy=E_energy+tmp_energy;
                
            end
        end
        
        tmp_num=sum(quad_w.*p_val(:,kx)*rho(Iu)*Jacobi_x/2);
        Num_part=Num_part+tmp_num;
        
    end
    
end



rho_max=max(rho(1:k:end));

E_momentum=E_momentum;
E_entropy=E_entropy;
E_energy=E_energy/2;

% [Num_part E_momentum E_entropy E_energy]
Eng=struct('num',Num_part,'mot',E_momentum,'entp',E_entropy,'eng',E_energy);

%***************************************
% Following is for Multiwavelet DG Basis
%***************************************

% Transfer to multi-DG bases
f_v=FMWT_COMP_v*f_v;
f_x=FMWT_COMP_x*f_x;


