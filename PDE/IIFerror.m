%This is for plotting the error between the error between exp(A*dt)*v and
%IIF exp

lev = 5; %level
deg = 4; %degree
interval = 5;
num = 20;%number of values of M tested
Mmax = interval*num; %maximum M

%dt = 2;
dt = 1/2^(lev); %time step
n = 5; %Number of time steps
%n = ceil(T/dt);
T = n*1/2^(lev);

%asgard(fokkerplanck1_4p4,'lev',lev,'deg',deg,'implicit',true,'num_steps',1)
%try diffusion2

%2d fokker planck
asgard(fokkerplanck2_6p1_withLHS,'lev',lev,'deg',deg,'implicit',true,'num_steps',1)
%asgard(fokkerplanck2_6p1,'lev',lev,'deg',deg,'implicit',true,'num_steps',1)

f0 = load('initial_vec000.mat','f0');
f0 = f0.f0; %f0 is initial vector

% load('pde.mat','pde')
% load('opts.mat','opts')
% load('hash_table.mat','hash_table')
% load('Meval.mat','Meval')
% load('nodes.mat','nodes')
load('matrix_iter000.mat','A');

AN = expm(A*T)*f0;
%an = wavelet_to_realspace(pde,opts,Meval,AN,hash_table);
M = linspace(interval,Mmax,num);
E = zeros(1,num);

[N,~] = size(A);
F = zeros(N,num);

for j=1:num
    m = interval*j;
    %IIF
    f = f0;
    for i=0:n-1
        [V,H] = myarnoldi(A,f,m);
        gamma = norm(f);
    
        kryExp = expm(H*dt);
        kryExp = kryExp(:,1); %first column in matrix
        kryExp = gamma*V*kryExp;

        f = kryExp;
    end
    %f1 = wavelet_to_realspace(pde,opts,Meval,f,hash_table); %Krylov 2nd Order
    F(:,j) = f;
    E(j) = norm(f-AN)/sqrt(N);
end

clf
semilogy(M,E);