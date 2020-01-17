x = linspace(100,500,401); %size of A
y = zeros(1,401); % error

%z = linspace(25,500,20); %M
texact = zeros(1,401);
tapprox = zeros(1,401);

M = 25;
%N = 500;

%A = rand(N);
%v = rand(N,1);
dt = 0.1;
%exact = expm(A*dt)*v;
I = eye(M);
e = I(:,1);

for i=1:401
A = rand(x(i));
v = rand(x(i),1);

tic
[V,H] = myarnoldi(A,v,M);
gamma = norm(v);
approx = gamma*V*expm(H*dt)*e;
tapprox(i) = toc;

tic
exact = expm(A*dt)*v;
texact(i) = toc;

%y(i) = norm(exact-approx); %error
end

%X = [ones(9,1),x'];

%y = log10(y);
%b = X\y';
%y2 = X*b;
tapprox = log10(tapprox);
texact = log10(texact);

%plot(x,y,x,y2)
%title('Logarithmic plot between error and size of A with M=25')
%xlabel('size of A') 
%ylabel('log of the error')
%legend({'error','slope ~0.023'},'Location','northwest')

plot(x,texact,x,tapprox)
title('Logarithmic plot between time and size of A with M=25')
xlabel('size of A') 
ylabel('log of the time')
legend({'exact','approximate'},'Location','northwest')