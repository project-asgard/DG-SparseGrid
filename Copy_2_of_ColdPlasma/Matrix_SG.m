function [GradMat,GradGradMat] = Matrix_SG(Lev,Deg,Lmax,pde)
%========================================================
% Construct the matrix for curl operator on [0,Lmax]
% Input:
%   Lev denotes the Level for the mesh
%   Deg denotes the degree for polynomial
%   Lmax denotes the domain for x variable
%   FMWT_COMP_x denotes the converting relationship
% Output:
%   GradMat:  \int u'v dx
%   GradGradMat: \int u'v' dx
%========================================================
close all
format short e


global A_encode %Mat Amat
%--DG parameters
quad_num=10;
%---------------

alpha = 10*Deg^3;
isGenAencode = 1;
isAssemble = 1;
c1 = 1;
c2 = 1;

% parameters for cg method
MaxIter = 1000;
Tol = 1e-6;


% compute the trace values
p_1 = legendre(-1,Deg);
p_2 = legendre( 1,Deg);

[quad_x,quad_w]=lgwt(quad_num,-1,1);
p_val = legendre(quad_x,Deg);
Dp_val = dlegendre(quad_x,Deg);
Dp_1 = dlegendre(-1,Deg);
Dp_2 = dlegendre(1,Deg);


%---------------------------
% Define Matrices
%---------------------------
nx = 2^(Lev);
hx = Lmax/nx;
dof_1D_x = Deg*nx;
dofs = 3*dof_1D_x^3;
GradMat = sparse(dof_1D_x,dof_1D_x);
GradGradMat = sparse(dof_1D_x,dof_1D_x);

%******************************************************
% generate 1D matrix for DG (du/dx,v) by weak form
% Here we only consider the central flux in the scheme
% -(u,dv/dx)+<{u},[v]>
%******************************************************
%======================================
% Matrices related to x variable GradX
%======================================
% compute (u',v)+1/2*u^{-}[v^{+}-v^{-}]
% generate 1D matrix for DG

% Matrix for (u',v')
val = Dp_val'*(quad_w*ones(1,Deg).*Dp_val)*1/hx/hx*2;
Ac = repmat({val},nx,1);
GG = blkdiag(Ac{:});

% Matrix for (u,v')
% take derivative for test function v
% val = Dp_val'*(quad_w*ones(1,Deg).*p_val)/hx;
% Ac = repmat({val},nx,1);
% G = blkdiag(Ac{:});


% Matrix for (u',v)
val = p_val'*(quad_w*ones(1,Deg).*Dp_val)/hx;
Ac = repmat({val},nx,1);
G = blkdiag(Ac{:});


% Matrix for (u,v) = eye


%****************************************
% Term for <{u_h},[v]>ds
%****************************************
Amd  = -p_1'*p_1/2+p_2'*p_2/2;
Asub = -p_1'*p_2/2;
Asup =  p_2'*p_1/2;
K = 1/hx*blktridiag([Amd],[Asub],[Asup],nx);
% Correction for boundary
K(1:Deg,1:Deg) = K(1:Deg,1:Deg)+(-p_1'*p_1)/hx/2;
K(end-Deg+[1:Deg],end-Deg+[1:Deg]) = K(end-Deg+[1:Deg],end-Deg+[1:Deg])+...
                (p_2'*p_2)/hx/2;

%****************************************
% Term for <[u_h],{v}>ds
%****************************************
Amd  = p_1'*(-p_1)/2+p_2'*p_2/2;
Asub =  p_1'*p_2/2;
Asup =  p_2'*(-p_1)/2;
H = 1/hx*blktridiag([Amd],[Asub],[Asup],nx);
% Correction for boundary
% H(1:Deg,1:Deg) = H(1:Deg,1:Deg)+(p_1'*(-p_1))/hx/2;
% H(end-Deg+[1:Deg],end-Deg+[1:Deg]) = H(end-Deg+[1:Deg],end-Deg+[1:Deg])+...
%                 (p_2'*p_2)/hx/2;

%****************************************
% Term for <{u_h'},[v]>ds
%****************************************
Amd = -p_1'*Dp_1/2+p_2'*Dp_2/2;
Asub =-p_1'*Dp_2/2;
Asup = p_2'*Dp_1/2;
L = 1/hx*blktridiag([Amd],[Asub],[Asup],nx)*2^(Lev+1);
% Correction for boundary
L(1:Deg,1:Deg) = L(1:Deg,1:Deg)+(-p_1'*Dp_1)/hx*2^(Lev+1)/2;
L(end-Deg+[1:Deg],end-Deg+[1:Deg]) = L(end-Deg+[1:Deg],end-Deg+[1:Deg])+...
                (p_2'*Dp_2)/hx*2^(Lev+1)/2;

%****************************************
% Term for <[u_h],{v'}>ds
%****************************************
Amd = Dp_1'*(-p_1)/2+Dp_2'*p_2/2;
Asub = Dp_1'*p_2/2;
Asup = Dp_2'*(-p_1)/2;
J = 1/hx*blktridiag([Amd],[Asub],[Asup],nx)*2^(Lev+1);
% Correction for boundary
% J(1:Deg,1:Deg) = J(1:Deg,1:Deg)+(Dp_1'*(-p_1))/hx*2^(Lev+1)/2;
% J(end-Deg+[1:Deg],end-Deg+[1:Deg]) = J(end-Deg+[1:Deg],end-Deg+[1:Deg])+...
%                 (Dp_2'*p_2)/hx*2^(Lev+1)/2;

%****************************************
% Term for <[u_h],[v]>ds
%****************************************
Amd  =  p_1'*(p_1)+p_2'*p_2;
Asub = -p_1'*p_2;
Asup =  p_2'*(-p_1);
Q = 1/hx*blktridiag([Amd],[Asub],[Asup],nx)*2/hx;

% convert these matries to the multiwavelet basis
FMWT_COMP = OperatorTwoScale(Deg,2^Lev);

GG = FMWT_COMP*GG*FMWT_COMP';
G = FMWT_COMP*G*FMWT_COMP';
H = FMWT_COMP*H*FMWT_COMP';
K = FMWT_COMP*K*FMWT_COMP';
L = FMWT_COMP*L*FMWT_COMP';
J = FMWT_COMP*J*FMWT_COMP';
Q = FMWT_COMP*Q*FMWT_COMP';

II = speye(dof_1D_x);

if isAssemble == 1
% tic
% For T1
t1 = kron(kron(II,GG),II);
t2 = kron(kron(II,II),GG);
t3 = kron(kron(GG,II),II);
tt = kron(G,G');
t4 = -kron(tt,II);
t5 = -kron(kron(G,II),G');
t6 = -kron(II,tt);

T1=[...
    t1+t2,t4,t5;...
    t4',t2+t3,t6;...
    t5',t6',t3+t1;...
    ];


% For T2

t1 = kron(kron(II,J),II);
t2 = kron(kron(II,II),J);
t3 = kron(kron(J,II),II);
tt = kron(H,G');
t4 = -kron(tt,II);
t5 = -kron(kron(H,II),G');
t6 = -kron(II,tt);
tt = kron(G',H);
t7 = -kron(tt,II);
t8 = -kron(kron(G',II),H);
t9 = -kron(II,tt);


T2=[...
    t1+t2,t4,t5;...
    t7,t2+t3,t6;...
    t8,t9,t3+t1;...
    ];


% For T3
t1 = kron(kron(II,L),II);
t2 = kron(kron(II,II),L);
t3 = kron(kron(L,II),II);
tt = kron(G,K);
t4 = -kron(tt,II);
t5 = -kron(kron(G,II),K);
t6 = -kron(II,tt);
tt = kron(K,G);
t7 = -kron(tt,II);
t8 = -kron(kron(K,II),G);
t9 = -kron(II,tt);


T3=[...
    t1+t2,t4,t5;...
    t7,t2+t3,t6;...
    t8,t9,t3+t1;...
    ];


% For T4
t1 = kron(kron(II,Q),II);
t2 = kron(kron(II,II),Q);
t3 = kron(kron(Q,II),II);


T4=blkdiag(t1+t2,t2+t3,t3+t1);
clear t1 t2 t3 t4 t5 t6 t7 t8 t9 tt

Mat = T1-c1*T2-c2*T3+alpha*T4-pde.w2*speye(dofs,dofs);
figure;subplot(2,2,1);spy(T1);subplot(2,2,2);spy(T2);subplot(2,2,3);spy(T3);subplot(2,2,4);spy(T4);
figure;spy(Mat);
111
end

A_encode = [];

if isGenAencode == 1
    % generate A_encode
    IndexI = [1:dof_1D_x^3];
    IndexJ = [1:dof_1D_x^3];
    
    % A11
    count = 1;
    A_encode{count}.IndexI = IndexI';
    A_encode{count}.IndexJ = IndexJ';
    
    A_encode{count}.A1=II;
    A_encode{count}.A2=GG-c1*J-c2*L+alpha*Q-pde.w2*II;
    A_encode{count}.A3=II;
%     full(sparse(IndexI',IndexJ',kron(II,kron(A_encode{count}.A2,II))))
    
    count = count+1;
    A_encode{count}.IndexI = IndexI';
    A_encode{count}.IndexJ = IndexJ';
    
    A_encode{count}.A1=II;
    A_encode{count}.A2=II;
    A_encode{count}.A3=GG-c1*J-c2*L+alpha*Q;
    
%     full(sparse(IndexI',IndexJ',kron(II,kron(II,A_encode{count}.A3))))
    
    % A12
    count = count+1;
    A_encode{count}.IndexI = IndexI';
    A_encode{count}.IndexJ = dof_1D_x^3+IndexJ';
    
    A_encode{count}.A1=-G;
    A_encode{count}.A2= G';
    A_encode{count}.A3=II;
    
    count = count+1;
    A_encode{count}.IndexI = IndexI';
    A_encode{count}.IndexJ = dof_1D_x^3+IndexJ';
    
    A_encode{count}.A1=H;
    A_encode{count}.A2=G';
    A_encode{count}.A3=c1*II;
    
    count = count+1;
    A_encode{count}.IndexI = IndexI';
    A_encode{count}.IndexJ = dof_1D_x^3+IndexJ';
    
    A_encode{count}.A1=G;
    A_encode{count}.A2=K;%H';
    A_encode{count}.A3=c2*II;
    
    % A13
    count = count+1;
    A_encode{count}.IndexI = IndexI';
    A_encode{count}.IndexJ = 2*dof_1D_x^3+IndexJ';
    
    A_encode{count}.A1=-G;
    A_encode{count}.A2=II;
    A_encode{count}.A3=G';
    
    count = count+1;
    A_encode{count}.IndexI = IndexI';
    A_encode{count}.IndexJ = 2*dof_1D_x^3+IndexJ';
    
    A_encode{count}.A1=H;
    A_encode{count}.A2=c1*II;
    A_encode{count}.A3=G';
    
    count = count+1;
    A_encode{count}.IndexI = IndexI';
    A_encode{count}.IndexJ = 2*dof_1D_x^3+IndexJ';
    
    A_encode{count}.A1=G;
    A_encode{count}.A2=c2*II;
    A_encode{count}.A3=K;%H';
    
    % A21
    count = count+1;
    A_encode{count}.IndexI = dof_1D_x^3+IndexI';
    A_encode{count}.IndexJ = IndexJ';
    
    A_encode{count}.A1=-G';
    A_encode{count}.A2=G;
    A_encode{count}.A3=II;
    
    count = count+1;
    A_encode{count}.IndexI = dof_1D_x^3+IndexI';
    A_encode{count}.IndexJ = IndexJ';
    
    A_encode{count}.A1=G';
    A_encode{count}.A2=H;
    A_encode{count}.A3=c1*II;
    
    count = count+1;
    A_encode{count}.IndexI = dof_1D_x^3+IndexI';
    A_encode{count}.IndexJ = IndexJ';
    
    A_encode{count}.A1=K;%H';
    A_encode{count}.A2=G;
    A_encode{count}.A3=c2*II;
    
    % A22
    count = count+1;
    A_encode{count}.IndexI = dof_1D_x^3+IndexI';
    A_encode{count}.IndexJ = dof_1D_x^3+IndexJ';
    
    A_encode{count}.A1=GG-c1*J-c2*L+alpha*Q-pde.w2*II;
    A_encode{count}.A2=II;
    A_encode{count}.A3=II;
    
    count = count+1;
    A_encode{count}.IndexI = dof_1D_x^3+IndexI';
    A_encode{count}.IndexJ = dof_1D_x^3+IndexJ';
    
    A_encode{count}.A1=II;
    A_encode{count}.A2=II;
    A_encode{count}.A3=GG-c1*J-c2*L+alpha*Q;
    
    % A23
    count = count+1;
    A_encode{count}.IndexI = dof_1D_x^3+IndexI';
    A_encode{count}.IndexJ = 2*dof_1D_x^3+IndexJ';
    
    A_encode{count}.A1=-II;
    A_encode{count}.A2=G;
    A_encode{count}.A3=G';
    
    count = count+1;
    A_encode{count}.IndexI = dof_1D_x^3+IndexI';
    A_encode{count}.IndexJ = 2*dof_1D_x^3+IndexJ';
    
    A_encode{count}.A1=c1*II;
    A_encode{count}.A2=H;
    A_encode{count}.A3=G';
    
    count = count+1;
    A_encode{count}.IndexI = dof_1D_x^3+IndexI';
    A_encode{count}.IndexJ = 2*dof_1D_x^3+IndexJ';
    
    A_encode{count}.A1=c2*II;
    A_encode{count}.A2=G;
    A_encode{count}.A3=K;%H';
    
    % A31
    count = count+1;
    A_encode{count}.IndexI = 2*dof_1D_x^3+IndexI';
    A_encode{count}.IndexJ = IndexJ';
    
    A_encode{count}.A1=-G';
    A_encode{count}.A2=II;
    A_encode{count}.A3=G;
    
    count = count+1;
    A_encode{count}.IndexI = 2*dof_1D_x^3+IndexI';
    A_encode{count}.IndexJ = IndexJ';
    
    A_encode{count}.A1=G';
    A_encode{count}.A2=c1*II;
    A_encode{count}.A3=H;
    
    count = count+1;
    A_encode{count}.IndexI = 2*dof_1D_x^3+IndexI';
    A_encode{count}.IndexJ = IndexJ';
    
    A_encode{count}.A1=K;%H';
    A_encode{count}.A2=c2*II;
    A_encode{count}.A3=G;
    
    % A32
    count = count+1;
    A_encode{count}.IndexI = 2*dof_1D_x^3+IndexI';
    A_encode{count}.IndexJ = dof_1D_x^3+IndexJ';
    
    A_encode{count}.A1=-II;
    A_encode{count}.A2=G';
    A_encode{count}.A3=G;
    
    count = count+1;
    A_encode{count}.IndexI = 2*dof_1D_x^3+IndexI';
    A_encode{count}.IndexJ = dof_1D_x^3+IndexJ';
    
    A_encode{count}.A1=c1*II;
    A_encode{count}.A2=G';
    A_encode{count}.A3=H;
    
    count = count+1;
    A_encode{count}.IndexI = 2*dof_1D_x^3+IndexI';
    A_encode{count}.IndexJ = dof_1D_x^3+IndexJ';
    
    A_encode{count}.A1=c2*II;
    A_encode{count}.A2=K;%H';
    A_encode{count}.A3=G;
    
    % A33
    count = count+1;
    A_encode{count}.IndexI = 2*dof_1D_x^3+IndexI';
    A_encode{count}.IndexJ = 2*dof_1D_x^3+IndexJ';
    
    A_encode{count}.A1=GG-c1*J-c2*L+alpha*Q-pde.w2*II;
    A_encode{count}.A2=II;
    A_encode{count}.A3=II;
    
    count = count+1;
    A_encode{count}.IndexI = 2*dof_1D_x^3+IndexI';
    A_encode{count}.IndexJ = 2*dof_1D_x^3+IndexJ';
    
    A_encode{count}.A1=II;
    A_encode{count}.A2=GG-c1*J-c2*L+alpha*Q;
    A_encode{count}.A3=II;
    
end



Lend = Lmax;Lstart = 0;
fx = zeros(quad_num,3);
fy = zeros(quad_num,3);
fz = zeros(quad_num,3);

hx = (Lend-Lstart)/nx;

for i=0:2^Lev-1
    
    xi_x = hx*(quad_x/2+1/2+i)+Lstart;
    fx = pde.rhs(xi_x,1/2+xi_x-xi_x,1/2+xi_x-xi_x);
    fy = pde.rhs(1/2+xi_x-xi_x,xi_x,1/2+xi_x-xi_x);
    fz = pde.rhs(1/2+xi_x-xi_x,1/2+xi_x-xi_x,xi_x);
    
    
    index = Deg*i+1:Deg*(i+1);
    f1x(index,1) = p_val'*(quad_w.*fx(:,1))*hx*sqrt(1/hx)/2;
    f1y(index,1) = p_val'*(quad_w.*fy(:,1))*hx*sqrt(1/hx)/2;
    f1z(index,1) = p_val'*(quad_w.*fz(:,1))*hx*sqrt(1/hx)/2;
    %     f1 = kron(kron(f1x,f1y),f1z);
    
    
    f2x(index,1) = p_val'*(quad_w.*fx(:,2))*hx*sqrt(1/hx)/2;
    f2y(index,1) = p_val'*(quad_w.*fy(:,2))*hx*sqrt(1/hx)/2;
    f2z(index,1) = p_val'*(quad_w.*fz(:,2))*hx*sqrt(1/hx)/2;
    %     f2 = kron(kron(f2x,f2y),f2z);
    
    f3x(index,1) = p_val'*(quad_w.*fx(:,3))*hx*sqrt(1/hx)/2;
    f3y(index,1) = p_val'*(quad_w.*fy(:,3))*hx*sqrt(1/hx)/2;
    f3z(index,1) = p_val'*(quad_w.*fz(:,3))*hx*sqrt(1/hx)/2;
    %     f3 = kron(kron(f3x,f3y),f3z);
    
end

f1x = FMWT_COMP*f1x;
f1y = FMWT_COMP*f1y;
f1z = FMWT_COMP*f1z;

f2x = FMWT_COMP*f2x;
f2y = FMWT_COMP*f2y;
f2z = FMWT_COMP*f2z;

f3x = FMWT_COMP*f3x;
f3y = FMWT_COMP*f3y;
f3z = FMWT_COMP*f3z;


ff =[kron(kron(f1x,f1y),f1z);kron(kron(f2x,f2y),f2z);kron(kron(f3x,f3y),f3z)]*(2*pi^2-pde.w2);

if isAssemble == 1

    sol = Mat\ff;


sol = Mat\ff;
% spy(Mat)
% full(Mat)

end

if isGenAencode == 1
x = zeros(dofs,1);

sol = cg(x,ff,MaxIter,Tol);
end


figure;plot(sol-ff/(2*pi^2-pde.w2))

uu = ff/(2*pi^2-pde.w2);
[max(abs(sol-uu)) norm(sol-uu)]
% compute the RHS term

% full([sol uu])
return


GradX = FMWT_COMP_x*GradX*FMWT_COMP_x';

end