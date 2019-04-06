% This code solves Helmholtz equation with DG-SG scheme

format short e

clc
clear
% close all
global PMat

addpath(genpath(pwd))
count = 1;
for Lev = 6:6
    
% Input for domain

Deg = 4;
Co = 1;

Ok = 3; Oh = 2; Op = 1;
% Deg = 2;
kappa = 1;%(Co*2^(Lev*Oh)*(Deg)^Op)^(1/Ok);
Lev = 9;

sigma =2;%1e-1i;
rho = 0;%1e-1i;
Ord = Deg-1; % Max is Deg-1

BCL = 'd';
BCR = 'r';


xMin = 0;
xMax = 1;
LMax = xMax - xMin;


% FunRHS = @(x)(5^2*sin(5*x)-kappa^2*sin(5*x));
% FunExt = @(x)(sin(5*x));
% FunBCL = @(x)(1i*kappa*sin(5*x)-5*cos(5*x));
% FunBCR = @(x)(1i*kappa*sin(5*x)+5*cos(5*x));

% FunRHS = @(x)(kappa^2*sin(kappa*x)-kappa^2*sin(kappa*x));
% FunExt = @(x)(sin(kappa*x));
% FunBCL = @(x)(1i*kappa*sin(kappa*x)-kappa*cos(kappa*x));
% FunBCR = @(x)(1i*kappa*sin(kappa*x)+kappa*cos(kappa*x));


% FunRHS = @(x)(x-x);
% FunExt = @(x)(exp(1i*kappa*x));
% FunGra = @(x)(1i*kappa*exp(1i*kappa*x));
% FunBCL = @(x)(1i*kappa*FunExt(x)-1*FunGra(x));
% FunBCR = @(x)(1i*kappa*FunExt(x)+1*FunGra(x));

FunRHS = @(x)(x-x-1);
FunExt = @(x)((1-cos(kappa*x)-sin(kappa)*sin(kappa*x)+1i*(1-cos(kappa))*sin(kappa*x))/kappa^2);
FunBCL = @(x)(x-x);
FunBCR = @(x)(x-x);

rhorho = 1;
% FunRHS = @(x)(kappa^2*sin(kappa*x)-rhorho*kappa^2*sin(kappa*x));
% FunExt = @(x)(sin(kappa*x));
% FunBCL = @(x)(x-x);
% FunBCR = @(x)(kappa*cos(kappa*x)+1i*kappa*sin(kappa*x));

N = 2^Lev; % # Grid Points
h = LMax/N;% size of the mesh
Jacobi = h/2;

% [Lev Deg kappa]% kappa^3*h^2/(Deg-1)^2]
%% Compute the 1D matrices
% compute the trace values for Polynomial and Derivatives
pL = legendre2(-1,Deg) * 1/sqrt(h);
pR = legendre2( 1,Deg) * 1/sqrt(h);
DDL = dlegendre2(-1,Deg,Deg-1) * 1/sqrt(h) * 2/h;
DDR = dlegendre2( 1,Deg,Deg-1) * 1/sqrt(h) * 2/h;
DDL = reshape(DDL,Deg,Deg)';
DDR = reshape(DDR,Deg,Deg)';

%---------------------------
% Matrices
%---------------------------
quad_num = Deg; % 10 Gaussian Quadratures

[quad_x,quad_w]=lgwt(quad_num,-1,1);
pVal  = legendre2(quad_x,Deg)  * 1/sqrt(h);
tmp = dlegendre2(quad_x,Deg,1);
DpVal = tmp(:,:,2) * 1/sqrt(h) * 2/h;


DoF_1D = Deg*N;
A = sparse(DoF_1D,DoF_1D);
b = sparse(DoF_1D,1);
bc = sparse(DoF_1D,1);
uu = sparse(DoF_1D,1);

% Volume Integral
VolInt = [DpVal'*(quad_w.*DpVal)]*Jacobi;

% Trace Integral
% B:= <{u'},[v]> + <[u],{v'}>
% [B1 B2 B3]
B1 = (-pL)' * (DDR(2,:)/2) + (DDL(2,:)'/2) * ( pR);
% B2L = (-pL)' * (DDL(2,:)/2)*2 + (DDL(2,:)'/2) * (-pL)*2+...
%     ( pR)' * (DDR(2,:)/2) + (DDR(2,:)'/2) * ( pR);
B2L =     ...(DDL(2,:)'/2) * (-pL)*2+...
    ( pR)' * (DDR(2,:)/2) + (DDR(2,:)'/2) * ( pR);
B2 = (-pL)' * (DDL(2,:)/2) + (DDL(2,:)'/2) * (-pL)+...
    ( pR)' * (DDR(2,:)/2) + (DDR(2,:)'/2) * ( pR);
% B2R = (-pL)' * (DDL(2,:)/2) + (DDL(2,:)'/2) * (-pL)+...
%       ( pR)' * (DDR(2,:)/2)*2+ (DDR(2,:)'/2) * ( pR)*2;
B2R = (-pL)' * (DDL(2,:)/2) + (DDL(2,:)'/2) * (-pL)+...
      0;%+ (DDR(2,:)'/2) * ( pR)*2;
B3 = ( pR)' * (DDL(2,:)/2) + (DDR(2,:)'/2) * (-pL);

% C:= I * p/h*<[u],[v]>
% [C1 C2 C3]
C1 = sigma*( (-pL)'*( pR) )* (Deg-1)/h;
C2L = sigma*( (pR)'*(pR) )* (Deg-1)/h;
C2 = sigma*( (-pL)'*(-pL)  + (pR)'*(pR) )* (Deg-1)/h;
C2R = sigma*( (-pL)'*(-pL) )* (Deg-1)/h;
C3 = sigma*( ( pR)'*(-pL) )* (Deg-1)/h;
% D:= I * sum_j (h/k)^(2*j-1)*<[D^ju],[D^jv]>
D1 = zeros(Deg,Deg);
D2 = zeros(Deg,Deg);
D3 = zeros(Deg,Deg);
D2L = zeros(Deg,Deg);
D2R = zeros(Deg,Deg);

for j = 1 : Ord%Deg-1
    coef = (h/(Deg-1))^(2*j-1);%*j;
    D1 = D1 + ( (-DDL(j+1,:))'*( DDR(j+1,:)) )*coef;
    D2 = D2 + ( (-DDL(j+1,:))'*(-DDL(j+1,:)) + DDR(j+1,:)'*DDR(j+1,:) )*coef;
    D2L = D2L + (   DDR(j+1,:)'*DDR(j+1,:) )*coef;
    D2R = D2R + ( (-DDL(j+1,:))'*(-DDL(j+1,:)) )*coef;
    D3 = D3 + ( ( DDR(j+1,:))'*(-DDL(j+1,:)) )*coef;
end


for Num = 0 : N-1
    xL = xMin + h*Num;
    xR = xL + h;
    
    % volume terms
    %     VolInt = VolInt;
    c = Deg*Num + [1:Deg];
    IndV = meshgrid(c);
    A = A + sparse(IndV',IndV,VolInt,DoF_1D,DoF_1D);
    
    if Num > 0 && Num < N-1
        % trace terms
        TraInt = [-B1+sigma*C1+rho*D1,...
                  -B2+sigma*C2+rho*D2,...
                  -B3+sigma*C3+rho*D3];
        
        IndU = [IndV-Deg,IndV,IndV+Deg];
        IndV = [IndV',IndV',IndV'];

    elseif Num == 0
%         TraInt = [-B2L+sigma*C2+rho*D2L, -B3+sigma*C3+rho*D3];
        IndU = [IndV, IndV+Deg];
        IndV = [IndV', IndV'];
        if BCL == 'r' % robin bc
            TraInt = [-B2L+sigma*C2L+rho*D2L + 1i*kappa*(-pL)'*(-pL), ...
                      -B3+sigma*C3+rho*D3];
            
            b(c) = b(c) + pL'*FunBCL(xMin);
        elseif BCL == 'd' %
            TraInt = [-B2L+sigma*C2L+rho*D2L-(-pL)' * (DDL(2,:)/2)*2-(DDL(2,:)'/2) * (-pL)*2, ...
                      -B3+sigma*C3+rho*D3];
            
            b(c) = b(c) + (DDL(2,:)')*FunExt(xMin);
        end

    elseif Num == N-1
%         TraInt = [-B1+sigma*C1+rho*D1, -B2R+sigma*C2+rho*D2R];
        
        IndU = [IndV-Deg, IndV];
        IndV = [IndV', IndV'];
        if BCR == 'r' % robin bc
            TraInt = [-B1+sigma*C1+rho*D1, ...
                      -B2R+sigma*C2R+rho*D2R + 1i*kappa*(pR)'*(pR)];
            
            b(c) = b(c) + pR'*FunBCR(xMax);
        elseif BCR == 'd' %
            TraInt = [-B1+sigma*C1+rho*D1, ...
            -B2R+sigma*C2R+rho*D2R-(pR)'*(DDR(2,:)/2)*2- (DDR(2,:)'/2) * ( pR)*2];
            
            b(c) = b(c) - (DDR(2,:)')*FunExt(xMax);
        end
    end
        A = A + sparse(IndU,IndV,TraInt,DoF_1D,DoF_1D);

    % assemble the Rhs
    xi = xMin + h *(quad_x/2+1/2+Num);
%     Val = pVal'*(quad_w.*sin(pi*xi))*Jacobi;
    Val = pVal'*(quad_w.*FunRHS(xi))*Jacobi;
    b(c) = b(c)+Val;
    uu(c) = pVal'*(quad_w.*FunExt(xi))*Jacobi;
end

M = A - rhorho*kappa^2*speye(DoF_1D,DoF_1D);

% M is the matrix

% CondM = condest(M);
% CondA = condest(A);


maxit = 1e5;

PMat = M-(1i/2)*kappa*speye(DoF_1D,DoF_1D);



% opts.my_shift =  epsilon;   % use variable opts to pass arguments to ?my_precond?
% gmres(A, b, restart, tol, maxit,  @my_precond, [],x0, opts)
 
x0 = b-b;
[sol,flag,relres,iter,resvec] = gmres(M,b,[],1e-6,maxit,PMat);

% MR = real(M);
% MI = imag(M);
% fR = real(b);
% fI =  imag(b);
% 
% MM = [MR  -MI; MI MR];
% bb = [fR;fI];
% % 
% [sol,info] = amg(MM,bb);
% 
% sol = sol(1:DoF_1D)+1i*sol(1+DoF_1D:end);
plot(resvec,'r-o')

% sol = M\b;


% v_M = eig(full(M));
% v_A = eig(full(A));
 
% figure;
% plot(v_M,'ro');hold on;plot(v_A,'b+');%plot(v_Ms,'k^');hold off;
% legend(['-\Delta-k^2 with Con = ',num2str(CondM)],...
%        ['-\Delta with Cond = ',num2str(CondA)]);
% %       ['ShiftMatrix with Cond= ',num2str(CondMs)])
% % % % % norm(sol-b)
% return
quad_num = Deg*3;
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


uh = Meval*(sol);
realuh = real(uh);
imaguh = imag(uh);

u = FunExt(x_node);
realu = real(u);
imagu = imag(u);

% [norm(uh-u) norm(sol-uu) sqrt(abs(sum(ww.*(uh-u).^2)))]
tmp = [Lev Deg kappa,kappa^Ok*h^Oh/(Deg-1)^Op...
    condest(M) ...
 ...sqrt(sum(ww.*(realuh-realu).^2+ww.*(imaguh-imagu).^2))]% ...
 sqrt(sum(ww.*(realuh-realu).^2+ww.*(imaguh-imagu).^2))/sqrt(sum(ww.*(realu).^2+ww.*(imagu).^2))]
r(count,:)=tmp;
count = count+1;
end
r

figure;
hold on
% subplot(1,2,1)
plot(x_node,realuh,'b-',x_node,imaguh,'r-','LineWidth',2)
% subplot(1,2,2)
% plot(x_node,realu,'b-',x_node,imagu,'r-','LineWidth',2) %-Meval*sol
% end