dt = 1e-5; %time step
lev = 5; %level
n = 15; %num_steps
deg = 4; %degree

asgard(mirror_velocity2,'timestep_method','BE', 'dt', dt, 'num_steps', 1, 'grid_type', 'SG', 'deg', deg, 'lev', lev)

M = 50;

f0 = load('initial_vec000.mat','f0');
f0 = f0.f0; %f0 is initial vector

load('pde.mat','pde')
load('opts.mat','opts')
load('hash_table.mat','hash_table')
load('Meval.mat','Meval')
load('nodes.mat','nodes')
load('matrix_iter000.mat','A');

[N,~] = size(A);
% e = eig(A);
% reale = real(e);
% m = max(reale)


f = f0;

for i=0:n-1
 
    [V,H] = myarnoldi(A,f,M);
    gamma = norm(f);
    
    kryExp = expm(H*dt);
    kryExp = kryExp(:,1); %first column in matrix
    kryExp = gamma*V*kryExp;

    f = kryExp;
end

g = f0;
I = eye(N);
for i = 0:n-1
    g = (I-dt*A)\g;
end

f0_real = wavelet_to_realspace(pde,opts,Meval,f0,hash_table);
F0 = reshape(f0_real,128,128);

an = expm(n*dt*A)*f0;
an_real = wavelet_to_realspace(pde,opts,Meval,an,hash_table);

f_real = wavelet_to_realspace(pde,opts,Meval,f,hash_table); %Krylov IIF

g_real = wavelet_to_realspace(pde,opts,Meval,g,hash_table); %Backward Euler

clf

plot_fval(pde,nodes,g_real,an_real,Meval) %analytic solution is calculated with matrix exponential.