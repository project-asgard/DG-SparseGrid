

clc
clear
% close all
format short e

kappa = 50;
Lev = 9;
N = 2^Lev;

xMin = 0;
xMax = 1;
LMax = xMax - xMin;

FunRHS = @(x)(x-x-1);
FunExt = @(x)((1-cos(kappa*x)-sin(kappa)*sin(kappa*x)+1i*(1-cos(kappa))*sin(kappa*x))/kappa^2);
FunBCL = @(x)(x-x);
FunBCR = @(x)(x-x);

N = 2^Lev; % # Grid Points
h = LMax/N;% size of the mesh
Jacobi = h/2;

A = sparse(N+1,N+1);
b = sparse(N+1,1);
Mass = sparse(N+1,N+1);
for j = 1 : N
    Mat_loc = 1/h*[1 -1;-1 1]-kappa^2*h*[1/3 1/6;1/6 1/3];
    Iu = [j,j,j+1,j+1];
    Iv = [j,j+1,j,j+1];
    A = A+sparse(Iu,Iv,Mat_loc(:),N+1,N+1);
    Mass = Mass + sparse(Iu,Iv,h*[1/3 1/6;1/6 1/3],N+1,N+1);
    
    x0 = xMin + h*(j-1);
    x1 = x0 + h;
    xi = (x0+x1)/2;
    
    Val = (FunRHS(xi))/2*h;
    b([j,j+1]) = b([j,j+1])+Val;
end
A(end,end) = A(end,end)+1i*kappa;
% A(end,

A(1,:) = 0;
A(1,1) = 1;
b(1) = FunExt(xMin);

% A(end,:) = 0;
% A(end,end) = 1;
% b(end) = 0;

% sol = A\b;


xx = [xMin:h:xMax];
% figure
% plot(xx,real(sol),'b--',xx,imag(sol),'r--','LineWidth',3)

M = A;
MR = real(M);
MI = imag(M);
fR = real(b);
fI =  imag(b);

MM = [MR  -MI; MI MR];
bb = [fR;fI];
maxit = 1e5;

PMat = M-(1i)/2*kappa^2*Mass;
PR = real(PMat);
PI = imag(PMat);

Pre = [PR  -PI; PI PR ];

% [x,flag,relres,iter,resvec] = bicg(MM,bb,1e-10,maxit,Pre);

% [x,flag,relres,iter,resvec] = pcg(MM,bb,1e-10,maxit,Pre);
% [x,flag,relres,iter,resvec] = gmres(MM,bb,[],1e-10,maxit,Pre);

% plot(resvec,'r-o')


[sol] = amg(Mat,ff);

[sol,flag,relres,iter,resvec] = gmres(M,b,[],1e-10,maxit,@MyFun);
figure
plot(xx,real(sol),'b-',xx,imag(sol),'r-','LineWidth',3)
% plot(resvec,'r-o')
