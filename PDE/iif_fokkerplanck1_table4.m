%for reproducing Table 4

lev = 5; %level
deg = 5; %degree


dt = 1/2^(lev);
%dt = 1/2^(2*lev); %time step
%dt = 5e-04;


asgard(fokkerplanck1_5p1a_noLHS,'lev',lev,'deg',deg,'implicit',true,'num_steps',1) %test 3.5.1a

M = 15; %krylov subspace approximation


f0 = load('initial_vec000.mat','f0');
f0 = f0.f0; %f0 is initial vector

load('pde.mat','pde')
load('opts.mat','opts')
load('hash_table.mat','hash_table')
load('Meval.mat','Meval')
load('nodes.mat','nodes')
load('matrix_iter000.mat','A');
[N,~] = size(A); %A is an N x N stiffness matrix


xx = nodes{1};

an = 4/sqrt(pi)*exp(-xx.^2); %Analytic solution

%%%% Krylov IIF%%%%
f = f0;
n_IIF=0;
f_prev = an;

tic
while norm(f-f_prev) > 1e-6
    f_prev = f;
    [V,H] = myarnoldi(A,f,M);
    gamma = norm(f);
    kryExp = expm(H*dt);
    kryExp = kryExp(:,1); %first column in matrix
    kryExp = gamma*V*kryExp;
    f = kryExp;
    n_IIF = n_IIF+1;
    %f1 = wavelet_to_realspace(pde,opts,Meval,g,hash_table);
    %plot(f1)
    %title(['time=',num2str(i)])
    %pause(0.0001)
end
toc
%%%%%%%%%%%%%%%%%%%
T_IIF = n_IIF*dt;

%%% Backwards Euler %%%%
h = f0;
I = eye(N);
n_BC = 0;
h_prev = an;
tic
while norm(h-h_prev) > 1e-6
    h_prev = h;
    h = (I-dt*A)\h;
    n_BC = n_BC+1;
end
toc
%%%%%%%%%%%%%%%%%%%%%%%%
T_BC = n_BC*dt;

f1 = wavelet_to_realspace(pde,opts,Meval,f,hash_table); %Krylov 2nd Order
h1 = wavelet_to_realspace(pde,opts,Meval,h,hash_table); %Backward Euler

%plot(xx,f1,'r-+',xx,an,'k.-','LineWidth',2) %K
plot(xx,f1,'r-h',xx,h1,'m-p',xx,an,'k.-','LineWidth',2) %K BE

%legend('Krylov IIF','Analytic Solution')
legend('Krylov IIF','Backward Euler','Analytic Solution')

title(['Fokker Planck Equation, T=',num2str(T_IIF),', dt=', num2str(dt,'%10.3e')])
norm(an-f1)
norm(an-h1,inf)
norm(an-h1)
%rms(an-f1)/rms(an)



%dt=1/2^lev, deg=1
%1.4142e-01 1.5222e-01
%4.6388e-02 5.2730e-02