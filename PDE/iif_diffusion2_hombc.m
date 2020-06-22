%for testing sourceless diffusion equation with homogeneous boundary condition.
%df/dt = d^2 f/dx^2 + d^2 f/dy^2 + s

T = 0.6; %end time
lev = 7; %level
deg = 3; %degree

%dt = 1/2^(lev); %time step
dt = 0.6;
%dt = 5e-06;

n = ceil(T/dt); %Number of time steps

asgard(diffusion2hombc,'lev',lev,'deg',deg,'implicit',true,'num_steps',1)

M = 5; %krylov subspace approximation

load('pde.mat','pde')
load('opts.mat','opts')
load('hash_table.mat','hash_table')
load('Meval.mat','Meval')
load('nodes.mat','nodes')

xx = nodes{1};
yy = nodes{2};
[XX,YY] = meshgrid(xx,yy);

%Analytic solution
an = exp(-T)*sin(2*pi*XX).*sin(2*pi*YY); %diffusion2

load('matrix_iter000.mat','A');
[N,~] = size(A); %A is an NxN stiffness matrix

f0 = load('initial_vec000.mat','f0');
f0 = f0.f0; %f0 is initial vector

%%%% Krylov IIF%%%%
f = f0;
tic
for i=0:n-1
    [V,H] = myarnoldi(A,f,M);
    gamma = norm(f);
    
    kryExp = expm(H*dt);
    kryExp = kryExp(:,1); %first column in matrix
    kryExp = gamma*V*kryExp;

    f = kryExp;
    
    %f1 = wavelet_to_realspace(pde,opts,Meval,g,hash_table);
    %plot(f1)
    %title(['time=',num2str(i)])
    %pause(0.0001)
end
toc
%%%%%%%%%%%%%%%%%%%
f1 = wavelet_to_realspace(pde,opts,Meval,f,hash_table); %Krylov 2nd Order
shape = sqrt(size(f1,1));

F = reshape(f1,shape,shape); %Krylov 2nd Order
format long e
rms(an(:)-F(:))

clf

contourf(F)
%%%% Forward Euler %%%%
g = f0;

%level 3, M=1 9.536290721172977e-02 M=5 9.768337318732324e-02 M=10 9.768337016472543e-02 M=100 9.768337016472545e-02
%level 4, 3.2506e-02
%level 5, 1.4851e-02
%level 6, 5.0056e-03
tic
for i=0:n-1   
    g = g + dt*A*g;
end
toc

%%% Backwards Euler %%%%
h = f0;
I = eye(N);
tic
for i=0:n-1
    h = (I-dt*A)\h;
end
toc
%%%%%%%%%%%%%%%%%%%%%%%%

g1 = wavelet_to_realspace(pde,opts,Meval,g,hash_table); %Forward Euler
h1 = wavelet_to_realspace(pde,opts,Meval,h,hash_table); %Backward Euler

%backwards euler
%level 3 9.473203082143739e-02
%level 5 1.407533141787036e-02

G = reshape(g1,shape,shape); %Forward Euler
H = reshape(h1,shape,shape); %Backward Euler

rms(an(:)-H(:))

surf(XX,YY,F)
%surf(XX,YY,P)
%surf(XX,YY,G)
%surf(XX,YY,H)

%surf(XX,YY,an)

%rms(an-p1)
%rms(an-f1)/rms(an)

title(['2D Diffusion Equation IIF, T=',num2str(T),', dt=', num2str(dt,'%10.3e')])