%for testing IIF method for 1D fokkerplanck equation

%T = 260; %end time
T = 2;
lev = 6; %level
deg = 4; %degree

%dt = 0.1;
dt = 1/2^(lev); %time step
%dt = 5e-05;

%0.177245

n = ceil(T/dt); %Number of time steps

%asgard(fokkerplanck1_5p1a_noLHS,'lev',lev,'deg',deg,'implicit',true,'num_steps',1) %test 3.5.1a
%asgard(fokkerplanck1_4p1a,'lev',lev,'deg',deg,'implicit',true,'num_steps',1) %test 3.2.3a
asgard(fokkerplanck1_4p4,'lev',lev,'deg',deg,'implicit',true,'num_steps',1) %Pitch angle collision
%asgard(fokkerplanck1_4p4,'lev',lev,'deg',deg,'implicit',true,'num_steps',1) %test 3.3

M = 141; %krylov subspace approximation

f0 = load('initial_vec000.mat','f0');
f0 = f0.f0; %f0 is initial vector

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
[N,~] = size(A); %A is an N x N stiffness matrix


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

%%%% Forward Euler %%%%
% g = f0;
% tic
% for i=0:n-1
%     g = g + dt*A*g;
%     
%     %%%for plotting %%%%
%     %g1 = wavelet_to_realspace(pde,opts,Meval,g,hash_table);
%     %plot(g1)
%     %title(['time=',num2str(i)])
%     %pause(0.001)
% end
% toc

%%% Backwards Euler %%%%
h = f0;
I = eye(N);
tic
for i=0:n-1
    h = (I-dt*A)\h;
end
toc
%%%%%%%%%%%%%%%%%%%%%%%%

xx = nodes{1};

f1 = wavelet_to_realspace(pde,opts,Meval,f,hash_table); %Krylov 2nd Order
%g1 = wavelet_to_realspace(pde,opts,Meval,g,hash_table); %Forward Euler
h1 = wavelet_to_realspace(pde,opts,Meval,h,hash_table); %Backward Euler



%Analytic/steady state solution

%phi = tanh(atanh(xx)-T);
%an = (1-phi.^2)./(1-xx.^2);
%an = 4/sqrt(pi)*exp(-xx.^2); %steady state solution to fokkerplanck1_5p1a
%A = 4;
%an = A / (2*sinh(A)) * exp(A*xx);
h = [3 0.5 1 0.7 3 0 3];
an = 0;
for L = 0:6
    an = an + h(L+1)*legendreP(L,xx)*exp(-L*(L+1)*T);%steady state solution to Fokkerplanck1_4p2
end


%plot(xx,f1,'r-+','LineWidth',2) %K
%plot(xx,f1,'r-+',xx,h1,'m-p','LineWidth',2) %K BE
%plot(xx,f1,'r-+',xx,an,'k.-','LineWidth',2) %K AN
%plot(xx,f1,'r-+',xx,g1,'b-o',xx,an,'k.-','LineWidth',2) %K FE AN
plot(xx,f1,'r-+',xx,h1,'m-p',xx,an,'k.-','LineWidth',2) %K BE AN
%plot(xx,f1,'r-+',xx,g1,'b-o',xx,h1,'m-p',xx,an,'k.-','LineWidth',2) %K FE BE AN

%legend('Krylov IIF')
%legend('Krylov IIF','Backward Euler')
%legend('Krylov IIF','Analytic/Steady State Solution')
%legend('Krylov IIF','Forward Euler','Analytic/Steady State Solution')
legend('Krylov IIF','Backward Euler','Analytic/Steady State Solution')
%legend('Krylov IIF','Forward Euler','Backwards Euler','Analytic/Steady State Solution')

title(['Fokker Planck Equation, T=',num2str(T),', dt=', num2str(dt,'%10.3e')])

rms(an-f1)
%rms(an-g1)
rms(an-h1)
%rms(an-f1)/rms(an)