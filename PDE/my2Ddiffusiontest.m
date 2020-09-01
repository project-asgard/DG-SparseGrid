%df/dt = d^2 f/dx^2 + d^2 f/dy^2 + s

T = 0.6; %end time
lev = 7; %level
deg = 2; %degree

dt = 1/2^(lev); %time step
%dt = 0.6;
%dt = 5e-06;

n = ceil(T/dt); %Number of time steps

asgard(diffusion2hombc,'lev',lev,'deg',deg,'implicit',true,'num_steps',1)

M = 100; %krylov subspace approximation

f0 = load('initial_vec000.mat','f0');
f0 = f0.f0; %f0 is initial vector

%bc0 = load('bc0.mat','bc0');
%bc0 = bc0.bc0;

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

shape = deg*2^lev;

F = reshape(f1,shape,shape); %Krylov 2nd Order
P = reshape(p1,shape,shape); %Krylov 3rd Order
G = reshape(g1,shape,shape); %Forward Euler
H = reshape(h1,shape,shape); %Backward Euler

xx = nodes{1};
yy = nodes{2};

[XX,YY] = meshgrid(xx,yy);

clf
surf(XX,YY,F)
%surf(XX,YY,P)
%surf(XX,YY,G)
%surf(XX,YY,H)

%Analytic solution
an = exp(-T)*sin(2*pi*XX).*sin(2*pi*YY); %diffusion2

%surf(XX,YY,an)

rms(an(:)-f1)
%rms(an-p1)
%rms(an-f1)/rms(an)

title(['2D Diffusion Equation IIF, T=',num2str(T),', dt=', num2str(dt,'%10.3e')])