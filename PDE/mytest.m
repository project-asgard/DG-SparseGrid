%asgard(fokkerplanck2_6p1_withLHS, 'implicit', true,'dt', 0.4, 'num_steps', 500, 'lev', 5, 'deg', 4, 'implicit_method', 'CN', 'time_independent_A', true)
%asgard(diffusion1, 'implicit', true,'dt', 0.4, 'num_steps', 500, 'lev', 5, 'deg', 4, 'implicit_method', 'CN', 'time_independent_A', true)
%asgard(advection1, 'implicit', true,'dt', 0.4, 'num_steps', 500, 'lev', 5, 'deg', 4, 'implicit_method', 'CN', 'time_independent_A', true)

n = 3000; %Number of time steps = ceil(T/dt)
dt = 2.5159e-02; %1.5e-2;

asgard(advection1,'lev',4,'deg',3,'implicit',true,'num_steps',1,'dt',dt)
%asgard(diffusion1,'lev',4,'deg',3,'implicit',true,'num_steps',n,'dt',dt)

M = 25; %krylov subspace approximation

filename2 = sprintf('initial_vec%03d.mat', 0);
f0 = load(filename2,'f0');
f0 = f0.f0; %f0 is initial vector

bc0 = load('bc0.mat','bc0');
bc0 = bc0.bc0;
s0 = load('solve_term','s0');
s = s0.s0;
s = s+bc0;

load('pde.mat','pde')
load('opts.mat','opts')
load('hash_table.mat','hash_table')
load('plotmatrix.mat','Meval','nodes')

%filename2 = sprintf('initial_vec%03d.mat', 1);
%f_2 = load(filename2,'f0');
%f_2 = f_2.f0;
    
filename1 = sprintf('matrix_iter%03d.mat', 0);
A = load(filename1,'A');
A = A.A;
[N,~] = size(A);

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
    
    f1 = wavelet_to_realspace(pde,opts,Meval,g,hash_table);
    plot(f1)
    title(['time=',num2str(i)])
    pause(0.001)
end
toc
%%%%%%%%%%%%%%%%%%%

%%%%% Krylov IIF 3rd Order %%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%I = eye(N);

%%%% Forward Euler %%%%
g = f0;
tic
for i=0:n-1
    g = g + dt*A*g+dt*s;
    
    g1 = wavelet_to_realspace(pde,opts,Meval,g,hash_table);
    plot(g1)
    title(['time=',num2str(i)])

    pause(0.001)
end
toc
%%%%%%%%%%%%%%%%%%%%%%%%


%%%% Backwards Euler %%%%
h = f0;
tic
for i=0:n-1
    h = (I-dt*A)\(h+dt*s);
end
toc
%%%%%%%%%%%%%%%%%%%%%%%%%


f1 = wavelet_to_realspace(pde,opts,Meval,f,hash_table);
g1 = wavelet_to_realspace(pde,opts,Meval,g,hash_table);
h1 = wavelet_to_realspace(pde,opts,Meval,h,hash_table);

xx = nodes{1};
an = cos(xx); %advection analytic solution
%an = cos(pi*x)*exp(-2*pi^2*n*dt); %diffusion analytic solution

lambda = eig(A);

% dtEstimate = (-lambda-conj(lambda))./(lambda.*conj(lambda));
% dtEstimate = dtEstimate(dtEstimate > 0);
% dtEstimate = min(dtEstimate);
%for lev = 4, deg = 3, dtEstimate = 2.5159e-02

%lambda = lambda*dt;
%plot(lambda,'ro','LineWidth',2)

%circle with radius 1 centered at -1
hold on
th = 0:pi/50:2*pi;
xunit = cos(th)-1;
yunit = sin(th);
plot(xunit, yunit);
%axis equal
hold off
%%%%%%%%%%%%%%%%%%%%%5

plot(xx,f1,'r-+',xx,g1,'b-o',xx,h1,'g-*',xx,an,'k.-','LineWidth',2)
%plot(xx,f1,'r-+',x,cos(x),'k.-','LineWidth',2)
%plot(xx,f1,'r-+',xx,g1,'b--','LineWidth',2)

%title('Krylov Advection, n=1000, dt=2.545e-02 lev=4, deg=3')
%xlabel('dt')
%ylabel('error')
legend('Krylov expm approximation','Forward Euler','Backwards Euler','f(x) = cos(x)')