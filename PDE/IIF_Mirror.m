dt = 1e-5; %time step
lev = 5; %level
n = 15; %num_steps
deg = 4; %degree

asgard(mirror_velocity2,'timestep_method','BE', 'dt', dt, 'num_steps', 15, 'grid_type', 'SG', 'deg', deg, 'lev', lev)

M = 50;

f0 = load('initial_vec000.mat','f0');
f0 = f0.f0; %f0 is initial vector

load('pde.mat','pde')
load('opts.mat','opts')
load('hash_table.mat','hash_table')
load('Meval.mat','Meval')
load('nodes.mat','nodes')
load('matrix_iter000.mat','A');
%load('bc0','bc0')

[N,~] = size(A);

f = f0;
% bc = zeros(N,n+1);
 %bc(:,1) = boundary_condition_vector(pde,opts,hash_table,0);
for i=0:n-1
    %bc(:,i+2) = boundary_condition_vector(pde,opts,hash_table,dt*(i+1));
    %[V,H] = myarnoldi(A,f+dt/2*bc(:,i+1),M);
    
    [V,H] = myarnoldi(A,f,M);
    gamma = norm(f);
    
    kryExp = expm(H*dt);
    kryExp = kryExp(:,1); %first column in matrix
    kryExp = gamma*V*kryExp;

    f = kryExp;
    %f = kryExp+dt/2*bc(:,i+2);
    
    
    %f1 = wavelet_to_realspace(pde,opts,Meval,g,hash_table);
    %plot(f1)
    %title(['time=',num2str(i)])
    %pause(0.0001)
end

g = f0;
I = eye(N);
for i = 0:n-1
    g = (I-dt*A)\g;
end

xx = nodes{1};
yy = nodes{2};
[XX,YY] = meshgrid(xx,yy);

f0_real = wavelet_to_realspace(pde,opts,Meval,f0,hash_table);
F0 = reshape(f0_real,128,128);

an = expm(n*dt*A)*f0;
an_real = wavelet_to_realspace(pde,opts,Meval,an,hash_table);
AN = reshape(an_real,128,128);

f_real = wavelet_to_realspace(pde,opts,Meval,f,hash_table); %Krylov IIF
F = reshape(f_real,128,128);

g_real = wavelet_to_realspace(pde,opts,Meval,g,hash_table); %Backward Euler

clf

%surf(XX,YY,AN)
%surf(XX,YY,F)
plot_fval(pde,nodes,an_real,g_real,Meval)
%contourf(F)
%contourf(AN)
