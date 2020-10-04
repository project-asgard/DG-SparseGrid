%df/dt = d^2 f/dx^2 - s
%soln_x = sin(pi*x);
%soln_t = exp(-2*t);

T = 0.1; %end time
lev = 4; %level
deg = 3; %degree

%dt = 1/2^(lev); %time step
dt = 5e-05;


n = ceil(T/dt); %Number of time steps

%asgard(diffusion1,'lev',lev,'deg',1,'implicit',true,'num_steps',1)
asgard(diffusion1,'lev',lev,'deg',deg,'implicit',true,'num_steps',1)
%asgard(fokkerplanck1_5p1a,'lev',4,'deg',3,'implicit',false,'num_steps',1)

%krylov subspace approximation
M = 25;


f0 = load('initial_vec000.mat','f0');
f0 = f0.f0; %f0 is initial vector

%bc0 = load('bc0.mat','bc0');
%bc0 = bc0.bc0;
%
%s0 = load('solve_term','s0');
%s = s0.s0;
%s = s+bc0;

load('pde.mat','pde')
load('opts.mat','opts')
load('hash_table.mat','hash_table')
load('Meval.mat','Meval')
load('nodes.mat','nodes')

%filename2 = sprintf('initial_vec%03d.mat', 1);
%f_2 = load(filename2,'f0');
%f_2 = f_2.f0;
    
%filename1 = sprintf('matrix_iter%03d.mat', 0);
load('matrix_iter000.mat','A');
%A = A.A;
[N,~] = size(A); %A is an NxN stiffness matrix

source0 = source_vector(pde,opts,hash_table,0);
bc0 = boundary_condition_vector(pde,opts,hash_table,0);

%%%% Krylov IIF %%%%
q = f0;
tic

for i=0:n-1
    [V,H] = myarnoldi(A,q,M);
    gamma = norm(q);
    
    kryExp = expm(H*dt);
    kryExp = kryExp(:,1); %first column in matrix
    kryExp = gamma*V*kryExp;

    q = kryExp;

    %f1 = wavelet_to_realspace(pde,opts,Meval,g,hash_table);
    %plot(f1)
    %title(['time=',num2str(i)])
    %pause(0.0001)
end

%%%% Krylov IIF 2nd order%%%%
f = f0;

tic
for i=0:n-1
    %source = source_vector(pde,opts,hash_table,dt*i);
    %bc = boundary_condition_vector(pde,opts,hash_table,dt*i);
    
    source = source0*exp(-2*dt*i);
    bc = bc0*exp(-2*dt*i);
    s = source+bc;
    f = f+dt/2*s;
    [V,H] = myarnoldi(A,f,M);
    gamma = norm(f);
    
    kryExp = expm(H*dt);
    kryExp = kryExp(:,1); %first column in matrix
    kryExp = gamma*V*kryExp;
    
    
    source2 = source0*exp(-2*dt*(i+1));
    bc2 = bc0*exp(-2*dt*(i+1));
    s2 = source2+bc2;
    f = kryExp+dt/2*s2;
    
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

    source = source0*exp(-2*dt*i);
    bc = bc0*exp(-2*dt*i);    
    s = source+bc;
    
    p_temp = p+(2/3)*dt*s;
    [V1,H1] = myarnoldi(A,p_temp,M);
    gamma1 = norm(p_temp);
    kryExp1 = expm(H1*dt);
    kryExp1 = kryExp1(:,1); %takes the first column
    kryExp1 = gamma1*V1*kryExp1;
    
    
    source3 = source0*exp(-2*dt*(i-1));
    bc3 = bc0*exp(-2*dt*(i-1));    
    s3 = source3+bc3;    
    
    [V2,H2] = myarnoldi(A,s3,M);
    gamma2 = norm(s3);
    kryExp2 = expm(2*H2*dt);
    kryExp2 = kryExp2(:,1); %first column
    kryExp2 = gamma2*V2*kryExp2;
    
    
    source2 = source0*exp(-2*dt*(i+1));
    bc2 = bc0*exp(-2*dt*(i+1));
    s2 = source2+bc2;

    p = kryExp1 - dt/12*kryExp2 + dt*(5/12)*s2;
end
toc

%%%% Forward Euler %%%%
g = f0;
tic
for i=0:n-1
    source = source0*exp(-2*dt*i);
    bc = bc0*exp(-2*dt*i);
    s = source+bc;
    
    g = g + dt*A*g+dt*s;
    
    %%%for plotting %%%%
    %g1 = wavelet_to_realspace(pde,opts,Meval,g,hash_table);
    %plot(g1)
    %title(['time=',num2str(i)])
    %pause(0.001)
end
toc

%%% Backwards Euler %%%%
h = f0;
I = eye(N);
tic
for i=0:n-1
    source2 = source0*exp(-2*dt*(i+1));
    bc2 = bc0*exp(-2*dt*(i+1));
    s2 = source2+bc2;

    h = (I-dt*A)\(h+dt*s2);
end
toc
%%%%%%%%%%%%%%%%%%%%%%%%

f1 = wavelet_to_realspace(pde,opts,Meval,f,hash_table); %Krylov 2nd Order
p1 = wavelet_to_realspace(pde,opts,Meval,p,hash_table); %Krylov 3rd Order
g1 = wavelet_to_realspace(pde,opts,Meval,g,hash_table); %Forward Euler
h1 = wavelet_to_realspace(pde,opts,Meval,h,hash_table); %Backward Euler
q1 = wavelet_to_realspace(pde,opts,Meval,q,hash_table); %Krylov IIF

xx = nodes{1};

%Analytic solution
%an = exp(-2*pi^2*dt*n)*cos(pi*xx); %diffusion1
an = exp(-T)*sin(xx); %diffusion1 alternative
%an = exp(-dt*n)*sin(2*pi*xx) ;%diffusionSourceless


%plot(xx,f1,'r-+',xx,an,'k.-','LineWidth',2) %K2
plot(xx,f1,'r-+',xx,q1,'m',xx,an,'k.-','LineWidth',2) %K2 K
%plot(xx,f1,'r-+',xx,g1,'b-o',xx,an,'k.-','LineWidth',2) %K2 FE
%plot(xx,f1,'r-+',xx,h1,'m-p',xx,an,'k.-','LineWidth',2) %K2 BE
%plot(xx,f1,'r-+',xx,p1,'m',xx,an,'k.-','LineWidth',2) %K2 K3
%plot(xx,f1,'r-+',xx,g1,'b-o',xx,h1,'m-p',xx,an,'k.-','LineWidth',2) %K2 FE BE
%plot(xx,f1,'r-+',xx,p1,'m',xx,h1,'g-*',xx,an,'k.-','LineWidth',2) %K2 K3 BE
%plot(xx,f1,'r-+',xx,p1,'m',xx,g1,'b-o',xx,h1,'g-*',xx,an,'k.-','LineWidth',2) %K2 K3 FE BE

%legend('Krylov IIF','Analytic Solution')
legend('Krylov IIF 2nd Order','Krylov IIF','Analytic Solution')
%legend('Krylov IIF','Forward Euler','Analytic Solution')
%legend('Krylov IIF','Backward Euler','Analytic Solution')
%legend('Krylov IIF 2nd order','Krylov IIF 3rd order','Analytic Solution')
%legend('Krylov IIF','Forward Euler','Backward Euler','Analytic Solution')
%legend('Krylov IIF 2nd order','Krylov IIF 3rd order','Backward Euler','Analytic Solution')
%legend('Krylov IIF 2nd order','Krylov IIF 3rd order','Forward Euler','Backwards Euler','Analytic Solution')

title(['Diffusion Equation, T=',num2str(T),', dt=', num2str(dt,'%10.3e')])

rms(an-f1)
%rms(an-p1)
%rms(an-f1)/rms(an)