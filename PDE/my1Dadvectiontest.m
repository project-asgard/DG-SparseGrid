%df/dt == -2*df/dx - 2*sin(x)

T = 1; %end time
%2.5159e-02
dt = 5e-3;

n = ceil(T/dt);

%asgard(advection1,'lev',4,'deg',3,'implicit',true,'num_steps',1,'dt',dt)
asgard(advection1,'lev',4,'deg',1,'implicit',true,'num_steps',1)

M = 2; %krylov subspace approximation

%filename2 = sprintf('initial_vec%03d.mat', 0);
f0 = load('initial_vec000.mat','f0');
f0 = f0.f0; %f0 is initial vector

bc0 = load('bc0.mat','bc0');
bc0 = bc0.bc0;
s0 = load('solve_term','s0');
s = s0.s0;
s = s+bc0;

%load('pde.mat','pde');
load('pde.mat','pde')
load('opts.mat','opts');
load('hash_table.mat','hash_table');
load('Meval.mat','Meval')
load('nodes.mat','nodes')

%filename2 = sprintf('initial_vec%03d.mat', 1);
%f_2 = load(filename2,'f0');
%f_2 = f_2.f0;
    
%filename1 = sprintf('matrix_iter%03d.mat', 0);
A = load('matrix_iter000.mat','A');
A = A.A;
[N,~] = size(A);

%%%% Forward Euler %%%%
g = f0;
tic
for i=0:n-1
    
    g = g + dt*A*g+dt*s;
    
    %%%% for plotting %%%%
    %g1 = wavelet_to_realspace(pde,opts,Meval,g,hash_table);
    %plot(g1)
    %title(['time=',num2str(i)])
    %pause(0.001)
end
toc
%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Krylov IIF 2nd order%%%%
f = f0;
tic
for i=0:n-1
    
    f = f+dt/2*s;
    [V,H] = myarnoldi(A,f,M);
    gamma = norm(f);
    %kryExp = gamma*V*expm(H*dt);
    %kryExp = kryExp(:,1); %first column in matrix
    kryExp = expm(H*dt);
    kryExp = kryExp(:,1); %first column in matrix
    kryExp = gamma*V*kryExp;
    
    f = kryExp+dt/2*s;
    
    %f1 = wavelet_to_realspace(pde,opts,Meval,g,hash_table);
    %plot(f1)
    %title(['time=',num2str(i)])
    %pause(0.0001)
end
toc
%%%%%%%%%%%%%%%%%%%

%%%%% Krylov IIF 3rd Order %%%%%%%
p = f0;
tic
for i=0:n-1
    
    p = p+(2/3)*dt*s;
    [V1,H1] = myarnoldi(A,p,M);
    gamma1 = norm(p);
    kryExp1 = expm(H1*dt);
    kryExp1 = kryExp1(:,1); %takes the first column
    kryExp1 = gamma*V1*kryExp1;
    
    [V2,H2] = myarnoldi(A,s,M);
    gamma2 = norm(s);
    kryExp2 = expm(2*H2*dt);
    kryExp2 = kryExp2(:,1); %first column
    kryExp2 = gamma2*V2*kryExp2;
        
    p = kryExp1 - dt/12*kryExp2 + dt*(5/12)*s;
 
%%%% for plotting %%%%    
%     p1 = wavelet_to_realspace(pde,opts,Meval,p,hash_table);
%     plot(xx,p1)
%     title(['time=',num2str(i)])
%     pause(0.001)
end
toc
%%%%%%%%%%%%%%%%%%%%%%%


%%%% Backwards Euler %%%%
h = f0;
I = eye(N);
tic
for i=0:n-1
    h = (I-dt*A)\(h+dt*s);
end
toc
%%%%%%%%%%%%%%%%%%%%%%%%%

f1 = wavelet_to_realspace(pde,opts,Meval,f,hash_table); %Krylov 2nd Order
g1 = wavelet_to_realspace(pde,opts,Meval,g,hash_table); %Forward Euler
h1 = wavelet_to_realspace(pde,opts,Meval,h,hash_table); %Backward Euler
p1 = wavelet_to_realspace(pde,opts,Meval,p,hash_table); %Krylov 3rd Order

xx = nodes{1};
an = cos(xx); %advection analytic solution

plot(xx,f1,'r-+',xx,an,'k.-','LineWidth',2) %K2
%plot(xx,f1,'r-+',xx,g1,'b-o',xx,an,'k.-','LineWidth',2) %K2 FE
%plot(xx,f1,'r-+',xx,h1,'m-p',xx,an,'k.-','LineWidth',2) %K2 BE
%plot(xx,f1,'r-+',xx,p1,'m',xx,an,'k.-','LineWidth',2) %K2 K3
%plot(xx,f1,'r-+',xx,g1,'b-o',xx,h1,'m-p',xx,an,'k.-','LineWidth',2) %K2 FE BE
%plot(xx,f1,'r-+',xx,p1,'m',xx,h1,'g-*',xx,an,'k.-','LineWidth',2) %K2 K3 BE
%plot(xx,f1,'r-+',xx,p1,'m',xx,g1,'b-o',xx,h1,'g-*',xx,an,'k.-','LineWidth',2) %K2 K3 FE BE

legend('Krylov IIF','Analytic Solution')
%legend('Krylov IIF','Forward Euler','Analytic Solution')
%legend('Krylov IIF','Backward Euler','Analytic Solution')
%legend('Krylov IIF 2nd order','Krylov IIF 3rd order','Analytic Solution')
%legend('Krylov IIF','Forward Euler','Backward Euler','Analytic Solution')
%legend('Krylov IIF 2nd order','Krylov IIF 3rd order','Backward Euler','Analytic Solution')
%legend('Krylov IIF 2nd order','Krylov IIF 3rd order','Forward Euler','Backwards Euler','Analytic Solution')

title(['Diffusion Equation, T=',num2str(T),', dt=', num2str(dt,'%10.3e')])