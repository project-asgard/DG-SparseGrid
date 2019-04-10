% This code solves 2D Helmholtz equation
% -Delta u - kappa^2 u = f
% u = g1, on Gamma_D
% gu/gn + i*kappa*u = g2, on Gamma_R
%========================================
clc
clear
close all
global A_encode
addpath(genpath(pwd))

Lev = 3;
Deg = 3;

DoF_1D = Deg * 2^Lev;
kappa = 1;

xMin = -pi; xMax = pi;
% FunRHS = @(x)(x-x-1);
% FunExt = @(x)((1-cos(kappa*x)-sin(kappa)*sin(kappa*x)+...
%                1i*(1-cos(kappa))*sin(kappa*x))/kappa^2);
% FunRHS = @(x)(pi^2*sin(pi*x));
% FunExt = @(x)(sin(pi*x));


% BC_opt.FunBCL = @(x)(x-x);
% BC_opt.FunBCR = @(x)(x-x);
% BC_opt.FunExt = FunExt;

FunRHS = @(x)(-kappa*sin(kappa*x));
FunExt = @(x)(sin(kappa*x));
BC_opt.FunBCL = @(x)(-kappa);
BC_opt.FunBCR = @(x)( kappa*cos(kappa*pi));
BC_opt.FunExt = FunExt;

% FunRHS = @(x)();
% FunExt = @(x)(sin(pi*x));
% BC_opt.FunBCL = @(x)(-kappa*cos(kappa*x)+1i*kappa*sin(kappa*x));
% BC_opt.FunBCR = @(x)( sin(kappa*x));
% BC_opt.FunExt = FunExt;

BC_opt.BCL = 'r';
BC_opt.BCR = 'd';

M = Helmholtz_1D(xMin,xMax,Lev,Deg,kappa,BC_opt);
% Mat = kron(M,speye(DoF_1D,DoF_1D))...
%     +kron(speye(DoF_1D,DoF_1D),M)...
%     -kappa^2*speye(DoF_1D^2,DoF_1D^2);

% % Creat HashTable
% Dim = 2;
% ComLev = perm_leq( Dim, Lev); % SG 
% ComLevIndex = [];



[b,bc] = RHSHelmholtz_1D(xMin,xMax,Lev,Deg,FunRHS,BC_opt);
 [coef] = ComputeCoef_1D(xMin,xMax,Lev,Deg,@(x)(sin(kappa*x)));

A_encode{1}.A = M;
A_encode{1}.B = speye(DoF_1D,DoF_1D);
A_encode{1}.IndexI = [1:DoF_1D^2];
A_encode{1}.IndexJ = [1:DoF_1D^2];

A_encode{2}.A = speye(DoF_1D,DoF_1D);
A_encode{2}.B = M;
A_encode{2}.IndexI = [1:DoF_1D^2];
A_encode{2}.IndexJ = [1:DoF_1D^2];

A_encode{3}.A = -kappa^2*speye(DoF_1D,DoF_1D);
A_encode{3}.B = speye(DoF_1D,DoF_1D);
A_encode{3}.IndexI = [1:DoF_1D^2];
A_encode{3}.IndexJ = [1:DoF_1D^2];

% A_encode is the matrix 

% The following is to assemble the global matrix
% AA = sparse(DoF_1D^2,DoF_1D^2);
% for i = 1:3
%     Mat_tmp = kron(A_encode{i}.A,A_encode{i}.B);
%     II = A_encode{i}.IndexI;
%     JJ = A_encode{i}.IndexJ;
%     [II,JJ] = meshgrid(II,JJ);
%     AA = AA + sparse(II,JJ,Mat_tmp,DoF_1D^2,DoF_1D^2);
% end



bb = kron(b,b)...
    -kron(bc,coef)-kron(coef,bc);


% B2 = ApplyA(bb);
maxit = 1e3;
% % 
[sol,flag,relres,iter,resvec] = gmres(@afun,bb,[],1e-3,maxit);
semilogy(resvec,'r-o')
% return
% sol_2D = AA\bb;

% % x1 = zeros(DoF_1D^2,1);
% % b1 = real(bb);
% % max_it = 1e8;
% % tol = 1e-10;
% % [sol_1,error,iter,flag]=cg_real(x1,b1,max_it,tol);
% % x1 = zeros(DoF_1D^2,1);
% % b1 = imag(bb);
% % '======= Imag Part ======'
% % [sol_2,error,iter,flag]=cg_imag(x1,b1,max_it,tol);
% % 
% % sol_cg = (sol_2+sol_1)/2+1i*(sol_2-sol_1)/2;

% sol = Mat\bb;

% plot
h = (xMax-xMin)/2^Lev;
quad_num = Deg;
[quad_x,quad_w]=lgwt(quad_num,-1,1);
ww = repmat(quad_w,2^Lev,1)*h/2;


Meval = sparse(quad_num*2^Lev,Deg*2^Lev);
x_node = zeros(quad_num*2^Lev,1);

for LL = 0:2^Lev-1
    Iu = Deg*LL+1:Deg*(LL+1); % row vector
    Iv = quad_num*LL+1:quad_num*(LL+1); % row vector
    
    Meval(Iv,Iu) = 1/sqrt(h)*legendre2(quad_x,Deg);
    %         Geval(Iv,Iu) = 2^(Lev/2)*dlegendre(quad_x,Deg)* 2*2^(Lev);
    
    xi = xMin+h*(quad_x/2+1/2+LL);
    x_node(Iv) = xi;
end

M2D = kron(Meval,Meval);

ureal = M2D*real(sol);uimag = M2D*imag(sol);
nz = size(x_node,1);

Fun2D = @(x,y)(sin(kappa*x).*sin(kappa*y));
[xx,yy] = meshgrid(x_node);

figure
subplot(2,2,1);
mesh(xx,yy,reshape(ureal,nz,nz),'FaceColor','interp','EdgeColor','none');
view(0,90);axis tight
subplot(2,2,2);
mesh(xx,yy,reshape(uimag,nz,nz),'FaceColor','interp','EdgeColor','none');
view(0,90);axis tight
% Fun2D = @(x,y)(...
%              ((1-cos(kappa*x)-sin(kappa)*sin(kappa*x)+...
%                1i*(1-cos(kappa))*sin(kappa*x)) ).*...
%              ((1-cos(kappa*y)-sin(kappa)*sin(kappa*y)+...
%                1i*(1-cos(kappa))*sin(kappa*y))/kappa^2));


u = Fun2D(xx,yy);
% figure;
subplot(2,2,3);
mesh(xx,yy,real(u),'FaceColor','interp','EdgeColor','none');
view(0,90);axis tight
subplot(2,2,4);
mesh(xx,yy,imag(u),'FaceColor','interp','EdgeColor','none');
view(0,90);axis tight
