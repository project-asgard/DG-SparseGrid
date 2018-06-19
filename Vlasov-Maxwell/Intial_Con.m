function OutPut=Intial_Con(Lev_x,Lev_v,k,Lmax,Vmax,pde,...
    FMWT_COMP_x,FMWT_COMP_v,Solver)
%==================================================================
% This code computes the initial conditions for f and rho=int_v fdv
% Input: Lev_x,Lev_v,k,Lmax,Vmax,pde,DimX,DimV
%   Here DimX is dimension for x, DimV is dimension for v
%       FMWT_COMP_x and FMWT_COMP_v are matrices for converting
% Output: f_v, f_x, rho, and Eng
%   Here f_v(:,DimV) and f_x(:,DimX)
%==================================================================
%--Quadrature
quad_num=10;
%---------------
[quad_x,quad_w]=lgwt(quad_num,-1,1);

p_val = legendre(quad_x,k);

%---------------------------
% Jacobi of variable x and v
%---------------------------
nx = 2^(Lev_x);hx = Lmax/nx;
Jacobi_x = hx;
dof_1D_x = k*nx;
f_x = zeros(dof_1D_x,pde.DimX);
if Solver == 'VM'
    E_0 = zeros(dof_1D_x,pde.E_Dim);
    B_0 = zeros(dof_1D_x,pde.B_Dim);
end


nv=2^(Lev_v);hv=2*Vmax/nv;
Jacobi_v=hv;
dof_1D_v=k*nv;
f_v=zeros(dof_1D_v,pde.DimV);

%=================================================
% Initial Condition for f_x
%=================================================
for Lx=0:nx-1
    xi_x=hx*(quad_x/2+1/2+Lx); % mapping from [-1,1] to physical domain
    %---------------------------------------------
    % Generate the coefficients for DG bases
    %---------------------------------------------
    Iu = [k*Lx+1:k*(Lx+1)];
    tmp_val = pde.Fx_0(xi_x);
    for dim = 1:pde.DimX
        f_x(Iu,dim)=p_val'*(quad_w.*tmp_val(:,dim))*Jacobi_x*sqrt(1/hx)/2;
    end
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
    tmp_val = pde.Fv_0(xi_v);
    for dim = 1:pde.DimV
        f_v(Iv,dim)=p_val'*(quad_w.*tmp_val(:,dim))*Jacobi_v*sqrt(1/hv)/2;
    end
end

%=================================================
% For Vlasov Maxwell Equation
%=================================================
if Solver == 'VM'
    
    for Lx=0:nx-1
        xi_x=hx*(quad_x/2+1/2+Lx); % mapping from [-1,1] to physical domain
        %---------------------------------------------
        % Generate the coefficients for DG bases
        %---------------------------------------------
        Iu = [k*Lx+1:k*(Lx+1)];
        tmp_Eval = pde.E_0(xi_x);
        tmp_Bval = pde.B_0(xi_x);
                
        for dim = 1:pde.E_Dim
            E_0(Iu,dim)=p_val'*(quad_w.*tmp_Eval(:,dim))*Jacobi_x*sqrt(1/hx)/2;
        end
        for dim = 1:pde.B_Dim
            B_0(Iu,dim)=p_val'*(quad_w.*tmp_Bval(:,dim))*Jacobi_x*sqrt(1/hx)/2;
        end
    end
end


%***************************************
% Following is for Multiwavelet DG Basis
%***************************************

% Transfer to multi-DG bases
f_v=FMWT_COMP_v*f_v;
f_x=FMWT_COMP_x*f_x;

if Solver == 'VM'
    E_0 = FMWT_COMP_x*E_0;
    B_0 = FMWT_COMP_x*B_0;
    OutPut = struct('f_v',f_v,'f_x',f_x,'E_0',E_0,'B_0',B_0);
else
    OutPut = struct('f_v',f_v,'f_x',f_x);
end

