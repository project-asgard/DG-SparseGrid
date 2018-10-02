% This code solves the 1D PDE with general B.C
% PDE is as follows::
%=============================================
%   d (au)/dx = f
% coupled with b.c
% u(0)=g(0) and/or u(1)=g(1)
% flux is chosen as
%   hat{f} = {{au}}+|a|*[[u]](1-\alpha)/2
% alpha (0 <= alpha <= 1)::
% 1 -- central flux
% 0 -- upwind flux
%=============================================

clc
clear
close all

format short e

% Test for generating 1D matrices
Lev = 5;
Deg = 2;
Np = 2^Lev;

dofs = Deg*Np;

Lmax = 1;
FuncA = @(x)(1);


BC_L  = 1:Deg;
BC_R = dofs-Deg+[1:Deg];


% RHS = @(x)(pi*cos(pi*x));
% exact = @(x)(sin(pi*x));
RHS = @(x)(-FuncA(0)*pi*sin(pi*x));
exact = @(x)(cos(pi*x));

alpha = 0;

%--Quadrature
quad_num=10;
%---------------



% compute the trace values
p_1 = legendre(-1,Deg);
p_2 = legendre(1,Deg);

[quad_x,quad_w]=lgwt(quad_num,-1,1);
p_val = legendre(quad_x,Deg);
Dp_val = dlegendre(quad_x,Deg);

nx=2^(Lev);
hx=Lmax/nx;


FluxMatrix = sparse(dofs,dofs);
MatrixStif = sparse(dofs,dofs);
b = sparse(dofs,1);

% basis with scaling functions
for num_cell = 0 : Np-1
    % cell = (xL,xR)
    xL = num_cell*hx; xR = xL+hx;
    
    % Term A*\int (u)(v')dx
    val_loc= FuncA((xL+xR)/2)*Dp_val'*(p_val.*quad_w)/hx;
    Index = ones(Deg,1)*[Deg*num_cell+1:Deg*(num_cell+1)];
    MatrixStif = MatrixStif+sparse(Index',Index,val_loc,dofs,dofs);
    
    xi = quad_x*(xR-xL)/2+(xR+xL)/2;
    val = sqrt(hx)/2*[p_val'*(quad_w.*RHS(xi))];
    b = b + sparse(Deg*num_cell+[1:Deg],ones(Deg,1),val,dofs,1);
    % Term (\hat{u})(v) on the surface
    % On the left end xL and right end xR
    
    vL = p_1; vR = p_2;
    uL_xL = p_2; uR_xL = p_1; % on xL
    uL_xR = p_2; uR_xR = p_1; % on xR
    nR = 1; nL = -1;
    aL = FuncA(xL); aR = FuncA(xR);
    c = Deg*num_cell+[1:Deg];
    
    if (num_cell >= 1) && (num_cell < Np-1)
        
        val_flux = [...
            nL*vL'*(aL+abs(aL)*(1-alpha)*nR)*uL_xL;...
            nL*vL'*(aL+abs(aL)*(1-alpha)*nL)*uR_xL;...
            nR*vR'*(aR+abs(aR)*(1-alpha)*nR)*uL_xR;...
            nR*vR'*(aR+abs(aR)*(1-alpha)*nL)*uR_xR;...
            ]*0.5/hx;
        IndexV = [c'*ones(1,Deg);c'*ones(1,Deg);c'*ones(1,Deg);c'*ones(1,Deg)];
        IndexU = [ones(Deg,1)*(c-Deg);ones(Deg,1)*(c);ones(Deg,1)*(c);ones(Deg,1)*(c+Deg) ];
        
    elseif num_cell == 0
        
        val_flux = [...
            nL*vL'*(aL+abs(aL)*(1-alpha)*nL)*uR_xL;...
            nR*vR'*(aR+abs(aR)*(1-alpha)*nR)*uL_xR;...
            nR*vR'*(aR+abs(aR)*(1-alpha)*nL)*uR_xR;...
            ]*0.5/hx;
        
        IndexV = [c'*ones(1,Deg);c'*ones(1,Deg);c'*ones(1,Deg)];
        IndexU = [ones(Deg,1)*(c);...
            ones(Deg,1)*(c);...
            ones(Deg,1)*(c+Deg) ];
        
        val_L = -nL*vL'*(aL+abs(aL)*(1-alpha)*nR)*exact(0)*0.5/sqrt(hx);
        
        
    elseif num_cell == Np-1
        
        val_flux = [...
            nL*vL'*(aL+abs(aL)*(1-alpha)*nR)*uL_xL;...
            nL*vL'*(aL+abs(aL)*(1-alpha)*nL)*uR_xL;...
            nR*vR'*(aR+abs(aR)*(1-alpha)*nR)*uL_xR;...
            ]*0.5/hx;
        
        IndexV = [c'*ones(1,Deg);c'*ones(1,Deg);c'*ones(1,Deg)];
        IndexU = [ones(Deg,1)*(c-Deg); ...
            ones(Deg,1)*(c);...
            ones(Deg,1)*(c)];
        
        val_R = -nR*vR'*(aR+abs(aR)*(1-alpha)*nL)*exact(1)*0.5/sqrt(hx);
        
    end
    
    FluxMatrix = FluxMatrix + sparse(IndexV,IndexU,val_flux,dofs,dofs);
    
    
end

Matrix = FluxMatrix - MatrixStif;

% Impose bc on the Input and Output side
b(BC_L) = b(BC_L) + val_L;
b(BC_R) = b(BC_R) + val_R;

x = Matrix\b;

num_plot = 1;
[quad_x,quad_w]=lgwt(num_plot,-1,1);
% quad_x = [-1,1]';

Lstart = 0;
Lend = 0;
p_val = legendre(quad_x,Deg);
for L=0:Np-1
    %---------------------------------------------
    % Generate the coefficients for DG bases
    %---------------------------------------------
    
    Iu = [Deg*L+1:Deg*(L+1)];
    
    Iv = [num_plot*L+1:num_plot*(L+1)];
    %     xi=h*(quad_x/2+1/2+L);
    
    x0 = Lstart+L*hx;
    x1 = x0+hx;
    xi = quad_x*(x1-x0)/2+(x1+x0)/2;[L*hx,L*hx+hx];%
    
    
    Meval(Iv,Iu)=sqrt(1/hx)*p_val;
    x_node(Iv,1)=xi;
    
end

% checked of projection
plot(x_node,Meval*x,'r-o',x_node,exact(x_node),'b--','LineWidth',2)

% figure;
% plot(x_node,Meval*b,'r-o',x_node,RHS(x_node),'b--');
