clc
clear
close all

format short e

% Test for generating 1D matrices
Lev = 1;
Deg = 2;
Np = 2^Lev;

%--Quadrature
quad_num=10;
%---------------

% compute the trace values
p_1 = legendre(-1,Deg);
p_2 = legendre(1,Deg);

[quad_x,quad_w]=lgwt(quad_num,-1,1);
p_val = legendre(quad_x,Deg);
Dp_val = dlegendre(quad_x,Deg);

load(['two_scale_rel_',num2str(Deg),'.mat'])
H1 = zeros(Deg);
G1 = zeros(Deg);

for j_x = 1:Deg
    for j_y = 1:Deg
        H1(j_x,j_y) = ((-1)^(j_x+j_y-2)  )*H0(j_x,j_y);
        G1(j_x,j_y) = ((-1)^(Deg+j_x+j_y-2))*G0(j_x,j_y);
    end
end

FMWT = zeros(Deg*Np);
iFMWT = zeros(Deg*Np);

k = Deg;
for j=1:Np/2
    % The reverse order from Lin
    FMWT(k*(j-1)+1:Deg*j,2*Deg*(j-1)+1:2*Deg*j)=[H0 H1];
    FMWT(k*(j+Np/2-1)+1:k*(j+Np/2),2*Deg*(j-1)+1:2*k*j) = [G0 G1];
end
iFMWT=FMWT';

sp = [];
FMWT_COMP = eye(k*Np);
n = Lev;
for j=1:Lev
    cFMWT = FMWT;
    % Reverse the index in matrix from Lin
    if j>1
        cFMWT = zeros(k*Np);
        cn = 2^(n-j+1)*k;
        cnr=Np*k-cn;
        cFMWT(cn+1:k*Np,cn+1:k*Np)=eye(Np*k-cn);
        cFMWT(1:cn/2,1:cn)=FMWT(1:cn/2,1:cn);
        cFMWT(cn/2+1:cn,1:cn)=FMWT(k*Np/2+1:k*Np/2+cn/2,1:cn);
    end

    FMWT_COMP = cFMWT*FMWT_COMP;
end

%---------------------------
% Jacobi of variable x 
% Define Matrices
%---------------------------
nx = 2^(Lev);
hx = 1/nx;
Jacobi_x = hx;
dof_1D_x = Deg*nx;


%-------------------------
% Tri-diagonal Matrix
%-------------------------
D =  Dp_val'*(quad_w.*p_val)...
    -p_1'*p_1/2+p_2'*p_2/2;% 1/hx
L = -p_1'*p_2/2; % 1/hx
U =  p_2'*p_1/2; % 1/hx

G = [G0 G1];
H = [H0 H1];

ZeroMat = zeros(Deg,Deg);

D1 = H*[D ZeroMat;ZeroMat D]*G';
L1 = H*[ZeroMat ZeroMat;L ZeroMat]*G';
U1 = H*[ZeroMat U;ZeroMat ZeroMat]*G';
Mat1 = D1+L1+U1;

D2 = G*[D ZeroMat;ZeroMat D]*H';
L2 = G*[ZeroMat ZeroMat;L ZeroMat]*H';
U2 = G*[ZeroMat U;ZeroMat ZeroMat]*H';
Mat2 = D2+L2+U2;

D3 = G*[D ZeroMat;ZeroMat D]*G';
L3 = G*[ZeroMat ZeroMat;L ZeroMat]*G';
U3 = G*[ZeroMat U;ZeroMat ZeroMat]*G';
Mat3 = D3+L3+U3;

A11 = D;
A12 = (Mat1)*2^1;
A13 = H*kron(Mat1,eye(2,2))
A21 = (Mat2)*2^1;
A22 = (Mat3)*2^1;

% % [H0*Mat1 H1*Mat1]

[A11 A12;A21 A22]

GradX = blktridiag([D/hx],[L/hx],[U/hx],nx);
% full(GradX)

GradX  = FMWT_COMP*GradX*FMWT_COMP'


for I = 1:Lev
    row = 2^(I-1)+1:2^I;
    G = blkdiag(G,G);
    H = blkdiag(H,H);
%     for J=1:Lev
%         col = 2^(J-1)+1:2^J;
%         tmp_A = 
%     end
end

figure;
subplot(1,2,1)
spy(G)
subplot(1,2,2)
spy(H)