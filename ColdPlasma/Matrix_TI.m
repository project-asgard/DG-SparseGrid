function [GradMat,GradGradMat] = Matrix_TI(Lev,Deg,Lmax,FMWT_COMP_x)
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
%--DG parameters
quad_num=10;
%---------------

alpha = 10;


% compute the trace values
p_1 = legendre(-1,Deg);
p_2 = legendre(1,Deg);

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
val = Dp_val'*(quad_w*ones(1,2).*Dp_val)*1/hx;
Ac = repmat({val},nx,1);
GG = blkdiag(Ac{:});

% Matrix for (u,v')
val = Dp_val'*(quad_w*ones(1,2).*p_val);
Ac = repmat({val},nx,1);
G = blkdiag(Ac{:});

% Matrix for (u',v)
val = p_val'*(quad_w*ones(1,2).*Dp_val);
Ac = repmat({val},nx,1);
G2 = blkdiag(Ac{:});


% Matrix for (u,v) = eye


%****************************************
% Term for <{u_h},[v]>ds
%****************************************
Amd  = -p_1'*p_1/2+p_2'*p_2/2;
Asub = -p_1'*p_2/2;
Asup =  p_2'*p_1/2;
K = 1/hx*blktridiag([Amd],[Asub],[Asup],nx);

%****************************************
% Term for <[u_h],{v}>ds
%****************************************
Amd  = p_1'*(-p_1)/2+p_2'*p_2/2;
Asub =  p_1'*p_2/2;
Asup =  p_2'*(-p_1)/2;
H = 1/hx*blktridiag([Amd],[Asub],[Asup],nx);

%****************************************
% Term for <{u_h'},[v]>ds
%****************************************
Amd = -p_1'*Dp_1/2+p_2'*Dp_2/2;
Asub =-p_1'*Dp_2/2;
Asup = p_2'*Dp_1/2;
L = blktridiag([Amd]/hx,[Asub]/hx,[Asup]/hx,nx);

%****************************************
% Term for <[u_h],{v'}>ds
%****************************************
Amd = Dp_1'*(-p_1)/2+Dp_2'*p_2/2;
Asub = Dp_1'*p_2/2;
Asup = Dp_2'*(-p_1)/2;
J = blktridiag([Amd]/hx,[Asub]/hx,[Asup]/hx,nx);

%****************************************
% Term for <[u_h],[v]>ds
%****************************************
Amd  =  p_1'*(p_1)+p_2'*p_2;
Asub = -p_1'*p_2/2;
Asup =  p_2'*(-p_1)/2;
Q = blktridiag([Amd]/hx,[Asub]/hx,[Asup]/hx,nx);

II = speye(dof_1D_x);
T1=[...
    kron(kron(II,GG),II)+kron(kron(II,II),GG),-kron(kron(G,G'),II),-kron(kron(G,II),G');...
    -kron(kron(G',G),II),kron(kron(GG,II),II)+kron(kron(II,II),GG),-kron(kron(II,G),G');...
    -kron(kron(G',II),G),-kron(kron(II,G'),G),kron(kron(GG,II),II)+kron(kron(II,GG),II);...
];
T2=[...
    kron(kron(II,J),II)+kron(kron(II,II),J),-kron(kron(H,G'),II),-kron(kron(H,II),G');...
    -kron(kron(G',H),II),kron(kron(J,II),II)+kron(kron(II,II),J),-kron(kron(II,H),G');...
    -kron(kron(G',II),H),-kron(kron(II,G'),H),kron(kron(J,II),II)+kron(kron(II,J),II);...
];
T3=[...
    kron(kron(II,L),II)+kron(kron(II,II),L),-kron(kron(G,K),II),-kron(kron(G,II),K);...
    -kron(kron(K,G),II),kron(kron(L,II),II)+kron(kron(II,II),L),-kron(kron(II,G),K);...
    -kron(kron(K,II),G),-kron(kron(II,K),G),kron(kron(L,II),II)+kron(kron(II,L),II);...
];
T4=[kron(kron(II,Q),II)+kron(kron(II,II),Q),zeros(dof_1D_x^3),zeros(dof_1D_x^3);...
    zeros(dof_1D_x^3),kron(kron(Q,II),II)+kron(kron(II,II),Q),zeros(dof_1D_x^3);...
    zeros(dof_1D_x^3),zeros(dof_1D_x^3),kron(kron(Q,II),II)+kron(kron(II,Q),II)];

Mat = T1-T2-T3+alpha*T4;
figure;subplot(2,2,1);spy(T1);subplot(2,2,2);spy(T2);subplot(2,2,3);spy(T3);subplot(2,2,4);spy(T4);
figure;spy(Mat);

% handling boundary condition

return


GradX = FMWT_COMP_x*GradX*FMWT_COMP_x';

end