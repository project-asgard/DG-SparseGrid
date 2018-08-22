% main code for LDG Poisson equation

% df/dt = d/dx(1-x^2)df/dx+source
% Method 1. LDG
% [A11 A12]
% [A21 A22]
% Here A11 = I, A12 = -(d/dx sqrt(1-x^2)q,p)
% A21 = -(sqrt(1-x^2)df/dx,w), A22 = (q,w) = I
% A12 = (sqrt(1-x^2)q,dp/dx)-<sqrt(1-x^2)\hat{q},p>
% A21 = (sqrt(1-x^2)f,dw/dx)-<sqrt(1-x^2)\hat{f},w>
% try with central flux
clear all
close all
clc

exactf = @(x,t)(exp(t)*sin(pi*x));
% source = @(x,t)(exp(t)*(sin(pi*x)+2*pi*x.*cos(pi*x)+(1-x.^2)*pi^2.*sin(pi*x)));
source = @(x)((sin(pi*x)+2*pi*x.*cos(pi*x)+(1-x.^2)*pi^2.*sin(pi*x)));
funcCoef = @(x)(sqrt(1-x.^2));

format short e
addpath(genpath(pwd))

Lev = 4;
Deg = 3;



Lstart = 0;
Lend = 1;
Lmax = Lend-Lstart;


%--Quadrature
quad_num=10;
%---------------

% compute the trace values
p_1 = legendre(-1,Deg);
p_2 = legendre(1,Deg);

[quad_x,quad_w]=lgwt(quad_num,-1,1);
p_val = legendre(quad_x,Deg);
Dp_val = dlegendre(quad_x,Deg);

%---------------------------
% Jacobi of variable x and v
% Define Matrices
%---------------------------
n=2^(Lev);h=Lmax/n;
Jacobi=h;
dof_1D=Deg*n;
A12 = sparse(dof_1D,dof_1D);
f0 = sparse(dof_1D,1);

b = sparse(dof_1D,1);

CFL = 0.01;
dt = CFL*h;


% generate 1D matrix for DG
for L=0:n-1

    %---------------------------------------------
    % (funcCoef*q,d/dx p)
    %---------------------------------------------
    x0 = Lstart+L*h;
    x1 = x0+h;
    xi = quad_x*(x1-x0)/2+(x1+x0)/2;
    
    val=1/h*[Dp_val'*(quad_w.*funcCoef(xi).*p_val)];
    c = Deg*L+1:Deg*(L+1);
    
    A12 = A12 + sparse(c'*ones(1,Deg),ones(Deg,1)*c,val,dof_1D,dof_1D);
    
    %----------------------------------------------
    % -<funcCoef*{q},p>
    %----------------------------------------------
    val=[ p_1'*funcCoef(x0)*p_2  p_1'*funcCoef(x0)*p_1,...
         -p_2'*funcCoef(x0)*p_2 -p_2'*funcCoef(x0)*p_1]/2;
    A12=A12+sparse(c'*ones(1,Deg),ones(Deg,1)*c,...
                   val(:,Deg+1:2*Deg)+val(:,2*Deg+1:3*Deg),...
                   dof_1D,dof_1D);

    if L>0
        A12=A12+sparse(c'*ones(1,Deg),ones(Deg,1)*c-Deg,val(:,1:Deg),dof_1D,dof_1D);
    end
    if L<n-1
        A12=A12+sparse(c'*ones(1,Deg),ones(Deg,1)*c+Deg,val(:,3*Deg+1:4*Deg),dof_1D,dof_1D);
    end
    
    val = sqrt(h)/2*[p_val'*(quad_w.*source(xi))]; 
    b(c)=val;
    
    val = sqrt(h)/2*[p_val'*(quad_w.*exactf(xi,0))]; 
    f0(c) = val;
    
    
end

% x = [Lstart:0.01:Lend];
% [f_loc] = EvalWavPoint4(Lstart,Lend,Lev,Deg,x,2);
% plot(f0)
[quad_x,quad_w]=lgwt(Deg,-1,1);
p_val = legendre(quad_x,Deg);
for L=0:n-1
    %---------------------------------------------
    % Generate the coefficients for DG bases
    %---------------------------------------------
    Iu=[Deg*L+1:Deg*(L+1)];
    xi=h*(quad_x/2+1/2+L);
    
    Meval(Iu,Iu)=sqrt(1/h)*p_val;
    x_node(Deg*L+1:Deg*L+Deg)=xi;

end

plot(x_node,Meval*f0)

% figure;plot(x,f_loc'*f0,'r-o');hold on;
% plot(x,exactf(x,time),'b--')

for t = 1:100
    time = t*dt;
    fval = f0+dt*A12*A12*f0+dt*b*exp(time);
    f0 = fval;
    
    plot(x_node,Meval*f0)
    pause
end
% figure;plot(x,f_loc'*f0,'r-o');hold on;
% plot(x,exactf(x,time),'b--')



