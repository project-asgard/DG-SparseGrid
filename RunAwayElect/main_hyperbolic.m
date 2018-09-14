% main code for LDG Poisson equation

% % df/dt = d/dx(1-x^2)df/dx+source
% df/dt = -d/dx ((1-x^2)f)
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
% clc

% Test
sigma = 0.1;
f0 = @(x)( exp(-x.^2/sigma^2) );
% f0 = @(x)(x-x+1);
phi = @(x,t)( tanh(atanh(x)-t) );
exactf = @(x,t)(...
    (1-phi(x,t).^2)./(1-x.^2).*f0(phi(x,t)) ...
    );
funcCoef = @(x)(1-x.^2);

format short e
addpath(genpath(pwd))


Lev = 6;
Deg = 2;
num_plot = Deg;
EndTime = 3;



Lstart = -1;
Lend = 1;
Lmax = Lend-Lstart;

FluxType = 'UF';

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

b = sparse(dof_1D,1);bb = sparse(dof_1D,1);
fexact = sparse(dof_1D,1);
qexact = sparse(dof_1D,1);

CFL = 0.01;
% dt = CFL*h^((Deg-1)/3)/2;
dt = CFL*h^((Deg)/3);

maxT = ceil(EndTime/dt)

% Assume
% [ I  A12]
% [A21  0 ]
% as the matrix
% e.g. we assume neuman boundary:: q=0 on boundary
%

A21 = sparse(dof_1D,dof_1D);
% generate 1D matrix for DG
for L=0:n-1
    
    %---------------------------------------------
    % (funcCoef*q,d/dx p)
    %---------------------------------------------
    x0 = Lstart+L*h;
    x1 = x0+h;
    xi = quad_x*(x1-x0)/2+(x1+x0)/2;
    xmid= (x0+x1)/2;
    
    val=1/h*[Dp_val'*(quad_w.*funcCoef(xi).*p_val)];
    
    c = Deg*L+1:Deg*(L+1);
    
    A12 = A12 + sparse(c'*ones(1,Deg),ones(Deg,1)*c,val,dof_1D,dof_1D);
    
    
    val = sqrt(h)/2*[p_val'*(quad_w.*exactf(xi,0))];
    fexact(c)=fexact(c)+val;
    
    
    %----------------------------------------------
    % -<funcCoef*{q},p>
    %----------------------------------------------
    if FluxType == 'CF'
        val=[ p_1'*funcCoef(x0)*p_2   p_1'*funcCoef(x0)*p_1,...
            -p_2'*funcCoef(x1)*p_2  -p_2'*funcCoef(x1)*p_1]/2/h;
        A12=A12+sparse(c'*ones(1,Deg),ones(Deg,1)*c,...
            val(:,Deg+1:2*Deg)+val(:,2*Deg+1:3*Deg),...
            dof_1D,dof_1D);
        
        
        if L>0
            A12=A12+sparse(c'*ones(1,Deg),ones(Deg,1)*c-Deg,val(:,1:Deg),dof_1D,dof_1D);
        end
        if L<n-1
            A12=A12+sparse(c'*ones(1,Deg),ones(Deg,1)*c+Deg,val(:,3*Deg+1:4*Deg),dof_1D,dof_1D);
        end
    end
    
    % up-winding flux
    if FluxType == 'UF'
        val=[p_1'*funcCoef(x0)*p_2   -p_2'*funcCoef(x1)*p_2]/h;
        A12=A12+sparse(c'*ones(1,Deg),ones(Deg,1)*c,...
            val(:,Deg+1:2*Deg),...
            dof_1D,dof_1D);
        
        
        if L>0
            A12=A12+sparse(c'*ones(1,Deg),ones(Deg,1)*c-Deg,val(:,1:Deg),dof_1D,dof_1D);
        end
    end
    
   if FluxType == 'DF'
        val=[p_1'*funcCoef(x0)*p_1   -p_2'*funcCoef(x1)*p_1]/h;
        A12=A12+sparse(c'*ones(1,Deg),ones(Deg,1)*c,...
            val(:,1:Deg),...
            dof_1D,dof_1D);
        
        
        if L<n-1
            A12=A12+sparse(c'*ones(1,Deg),ones(Deg,1)*c+Deg,val(:,1+Deg:2*Deg),dof_1D,dof_1D);
        end
    end
    
    
    
    val = sqrt(h)/2*[p_val'*(quad_w.*exactf(xi,0))];
    f0(c) = val;
    
    
end



[quad_x,quad_w]=lgwt(num_plot,-1,1);
% quad_x = [-1,1]';

p_val = legendre(quad_x,Deg);
for L=0:n-1
    %---------------------------------------------
    % Generate the coefficients for DG bases
    %---------------------------------------------
    
    Iu = [Deg*L+1:Deg*(L+1)];
    
    Iv = [num_plot*L+1:num_plot*(L+1)];
    %     xi=h*(quad_x/2+1/2+L);
    
    x0 = Lstart+L*h;
    x1 = x0+h;
    xi = quad_x*(x1-x0)/2+(x1+x0)/2;[L*h,L*h+h];%
    
    
    Meval(Iv,Iu)=sqrt(1/h)*p_val;
    x_node(Iv,1)=xi;
    
end

% checked of projection
plot(x_node,Meval*f0,'r-o',x_node,exactf(x_node,0),'b--','LineWidth',2)

% return
Mat = A12;
% return
% max(abs(Meval*f0))
% val = Meval*A12*f0-exactq(x_node,0);
% [norm(val) max(abs(val))]

b = b+bb;


[quad_x,quad_w]=lgwt(num_plot,-1,1);
total_particle = 0;
ffval = Meval*f0;
L2_stability = 0;
for i = 1:num_plot
    total_particle =  total_particle+quad_w(i)*h/2*sum(ffval(i:num_plot:end));
    L2_stability = L2_stability+quad_w(i)*h/2*sum(ffval(i:num_plot:end).^2);
end
[total_particle L2_stability]

tp(1) = total_particle;
Lp(1) = L2_stability;

figure
for t = 1:maxT
    %     t
    time = t*dt;
    %     tmp = A12'*A12*f0;%dt*A12'*A12*f0+dt*b*exp(time);
    %     tmp2 = A12*A12*f0;
    
    % %     fval = f0+dt*Mat*f0+dt*b;%*exp(time);
    
    f1 = f0 + dt*( Mat*f0+b*exp(time-dt) );
    f2 = 3/4*f0+1/4*f1+1/4*dt*(Mat*f1+b*exp(time));
    fval = 1/3*f0+2/3*f2+2/3*dt*(Mat*f2+b*exp(time-dt/2));
    
    %     f1 = f0 + dt*( Mat*f0+b );
    %     f2 = 3/4*f0+1/4*f1+1/4*dt*(Mat*f1+b );
    %     fval = 1/3*f0+2/3*f2+2/3*dt*(Mat*f2+b );
    
    f0 = fval;
    
    plot(x_node,Meval*f0,'r-o',x_node,Meval*(A12*f0),'b-<')
    title(['time at ',num2str(time)])
    pause (0.1)
    
    total_particle = 0;
    L2_stability = 0;
    ffval = Meval*f0;
    for i = 1:num_plot
        
        total_particle =  total_particle+...
            quad_w(i)*h/2*sum(ffval(i:num_plot:end));
        L2_stability = L2_stability+quad_w(i)*h/2*sum(ffval(i:num_plot:end).^2);
    end
    tp(t+1) = total_particle;
    Lp(t+1) = L2_stability;
    
    if abs(time-0.5)<=dt || abs(time-1)<=dt || abs(time-2)<=dt || abs(time-3)<=dt
        save(['hyper_',FluxType,'_Deg',num2str(Deg),'_Lev',num2str(Lev),'_End',num2str(time),'.mat'])
    end
end
% figure;plot(x,f_loc'*f0,'r-o');hold on;
% plot(x,exactf(x,time),'b--')
% hold on
% plot(x_node,exactf(x_node,time),'r-o')
%  max(abs(Meval*f0))
val = Meval*f0-exactf(x_node,time);
%  [norm(val) max(abs(val))]

fL2 = 0; fLinf = max(abs(val));
%  total_particle = 0;
ffval = Meval*f0;
for i = 1:num_plot
    fL2 = fL2 + quad_w(i)*h/2*sum(val(i:num_plot:end).^2);
    %       total_particle =  total_particle+quad_w(i)*h/2*sum(ffval(i:num_plot:end));
end
[sqrt(fL2) fLinf]
%  total_particle


%  err = f0-fexact;
%  full([norm(err) max(abs(err))])
%
%  err = A12*f0-qexact;
%  full([norm(err) max(abs(err))])

figure;
plot(x_node,Meval*f0,'r-o',x_node,exactf(x_node,time),'r--','LineWidth',2);


figure;
plot(tp,'r-o'); hold on; plot(Lp,'b-o')
