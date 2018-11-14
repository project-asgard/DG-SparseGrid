% function [GradMat,GradGradMat] = Matrix_SG2(Lev,Deg,Lmax,pde)
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
clear 
<<<<<<< HEAD
Lev = 2;
=======
Lev = 5;
>>>>>>> 195f27eea2d0548bb1db5ce8dac04acc7219e171
Deg = 2;
Lmax = 1;
pde = Maxwell1;

AssembleMatrix = 0;
% Lev = 1;Deg = 1;Lmax = 1;pde=Maxwell1;

global A_encode invM%Mat Amat
%--DG parameters
quad_num=10;
%---------------

<<<<<<< HEAD
alpha = 2*(Deg + 1)*(Deg + 3);%1000*Deg^3;
=======
alpha = 10*Deg*(Deg+1);
>>>>>>> 195f27eea2d0548bb1db5ce8dac04acc7219e171
isGenAencode = 1;
isAssemble = 0;
c1 = 1;
c2 = 1;

% parameters for cg method
MaxIter = 3000;
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
%%
H(1:Deg,1:Deg) = H(1:Deg,1:Deg)+(p_1'*(-p_1))/hx/2;
H(end-Deg+[1:Deg],end-Deg+[1:Deg]) = H(end-Deg+[1:Deg],end-Deg+[1:Deg])+...
                (p_2'*p_2)/hx/2;

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
%%
J(1:Deg,1:Deg) = J(1:Deg,1:Deg)+(Dp_1'*(-p_1))/hx*2^(Lev+1)/2;
J(end-Deg+[1:Deg],end-Deg+[1:Deg]) = J(end-Deg+[1:Deg],end-Deg+[1:Deg])+...
                (Dp_2'*p_2)/hx*2^(Lev+1)/2;

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

% Generate Hash Table
[forwardHash,inverseHash]=HashTable(Lev,3);

% Assemble the global matrix
Dofs = forwardHash.dof;
Dofs3 = 3*Dofs;

combs = perm_leq(3,Lev);

% combs = [0,0,0;0 0 1;0 1 0;0 1 1;1 0 0;1 0 1;1 1 0;1 1 1];
% Dofs = size(combs,1);
% Dofs3 = 3*Dofs;

% full grid
val = [0:Lev];
n3 = repmat(val,1,length(val)*length(val));
n3 = n3(:);
n2 = repmat(val,length(val),length(val));
n2 = n2(:);
n1 = repmat(val,length(val)*length(val),1);
n1 = n1(:);

combs = [n1,n2,n3];

combs_index = zeros(size(combs,1),2);
count = 0;
for l = 1:size(combs,1)
    n1 = combs(l,1);
    n2 = combs(l,2);
    n3 = combs(l,3);
    
    nz1 = 2^max(n1-1,0);
    nz2 = 2^max(n2-1,0);
    nz3 = 2^max(n3-1,0);
    
    count_loc = nz1*nz2*nz3;
    combs_index(l,:) = [count+1,count+count_loc];
    count = count+count_loc;
end
Dofs = count;
Dofs3 = 3*Dofs;
% [combs combs_index]

% clear A_encode
A_encode=[];

count= 0;
% generate the combination of lev
for l1 = 1:size(combs,1)

    row1 = Deg*(2^combs(l1,1)-2^max(0,combs(l1,1)-1))+1:Deg*2^combs(l1,1);
    row2 = Deg*(2^combs(l1,2)-2^max(0,combs(l1,2)-1))+1:Deg*2^combs(l1,2);
    row3 = Deg*(2^combs(l1,3)-2^max(0,combs(l1,3)-1))+1:Deg*2^combs(l1,3);
    
    IndexI = [Deg^3*(combs_index(l1,1)-1)+1:Deg^3*combs_index(l1,2)];

%     tmp = permute(reshape(IndexI,),[2,1])
    
    for l2 = 1:size(combs,1)

        col1 = Deg*(2^combs(l2,1)-2^max(0,combs(l2,1)-1))+1:Deg*2^combs(l2,1);
        col2 = Deg*(2^combs(l2,2)-2^max(0,combs(l2,2)-1))+1:Deg*2^combs(l2,2);
        col3 = Deg*(2^combs(l2,3)-2^max(0,combs(l2,3)-1))+1:Deg*2^combs(l2,3);

        
        IndexJ = [Deg^3*(combs_index(l2,1)-1)+1:Deg^3*combs_index(l2,2)];
        
%         IndexI
%         IndexJ
        A_loc = AssembleAencodeMaxwell(IndexI,IndexJ,GG,G,J,H,K,L,Q,pde.w2,Deg^3*Dofs,row1,row2,row3,col1,col2,col3,alpha);
%         Mat = Aencode2Matrix(A_loc,Deg^3*Dofs3);
%         spy(Mat(1:16,1:16))

        if size(A_loc,2)>0 % only add the non-zero A_encode
            for nz = 1:size(A_loc,2)
                A_encode{count+nz}.IndexI = A_loc{nz}.IndexI;
                A_encode{count+nz}.IndexJ = A_loc{nz}.IndexJ;
                A_encode{count+nz}.A1 = A_loc{nz}.A1;
                A_encode{count+nz}.A2 = A_loc{nz}.A2;
                A_encode{count+nz}.A3 = A_loc{nz}.A3;
  
            end
        end
        count = count+size(A_loc,2);

    end
    
    
end

size(A_encode)

% Mat = Aencode2Matrix(A_encode,Deg^3*Dofs3);
% % figure;
% % tmp = Mat-Mat';
% % tmp(find(abs(tmp)<1e-5))=0;
% % % spy(Mat)
% % spy(tmp)

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

% Assemble the RHS term for Sparse Grids
ff = sparse(Deg^3*Dofs3,1);
for l1 = 1:size(combs,1)

    row1 = Deg*(2^combs(l1,1)-2^max(0,combs(l1,1)-1))+1:Deg*2^combs(l1,1);
    row2 = Deg*(2^combs(l1,2)-2^max(0,combs(l1,2)-1))+1:Deg*2^combs(l1,2);
    row3 = Deg*(2^combs(l1,3)-2^max(0,combs(l1,3)-1))+1:Deg*2^combs(l1,3);
    
    IndexI = [Deg^3*(combs_index(l1,1)-1)+1:Deg^3*combs_index(l1,2)];
%     [row1 row2 row3]
    
    ff(IndexI) = kron(kron(f1x(row1),f1y(row2)),f1z(row3));
    ff(Deg^3*Dofs+IndexI) = kron(kron(f2x(row1),f2y(row2)),f2z(row3));
    ff(2*Deg^3*Dofs+IndexI) = kron(kron(f3x(row1),f3y(row2)),f3z(row3));
end


% ff =[kron(kron(f1x,f1y),f1z);kron(kron(f2x,f2y),f2z);kron(kron(f3x,f3y),f3z)]*(2*pi^2-pde.w2);
ff = ff*(2*pi^2-pde.w2);



if AssembleMatrix == 1
Mat = Aencode2Matrix(A_encode,Deg^3*Dofs3);
% tmp = Mat-Mat';
% tmp(find(abs(tmp)<1e-5))=0;
% spy(tmp)
% condest(Mat)
% figure;
% spy(Mat)
% % return
sol = Mat\ff;
uu = ff/(2*pi^2-pde.w2);
% % full([max(abs(sol-uu)) norm(sol-uu)/norm(uu)])
full([max(abs(sol-uu)) norm(sol-uu)])
% figure
% plot(sol-uu)
% return
end

% [ 'Mat_Lap.m'];

% Lap = GG-L-J+alpha*Q;
% 
% dofs = size(Lap,1);
% II = speye(dofs);
% 
% Mat_Lap = kron(II,kron(II,Lap))+kron(II,kron(Lap,II))+kron(Lap,kron(II,II));
% 
% invMat = inv(Mat_Lap);
% 
% invM = blkdiag(invMat,invMat,invMat);

% % Mat = Aencode2Matrix(A_encode,Deg^3*Dofs3);
% % MM = diag(diag(Mat));
% % invM = inv(MM);



x = zeros(Deg^3*Dofs3,1);%
sol = cg(x,ff,MaxIter,Tol);

% Mat = Aencode2Matrix(A_encode,Deg^3*Dofs3);
% %  [sol] = amg(Mat,ff);
% sol3 = Mat\ff;
% [sol,flag2,rr2,iter2,rv2] = pcg(@afun,ff,Tol,MaxIter);

% % plot(sol2-sol)

% tol = 1e-6;  maxit = 1000; 
% sol3 = gmres(@afun,ff,10,Tol,MaxIter);
% % 
% [sol3,fl0,rr0,it0,rv0]= gmres(@afun,ff,10,Tol,MaxIter);
% semilogy(0:length(rv0)-1,rv0/norm(ff),'-o');
% xlabel('Iteration number');
% ylabel('Relative residual');
% 
% % % figure;plot(sol2-sol)



uu = ff/(2*pi^2-pde.w2);
[max(abs(sol-uu)) norm(sol-uu)]
% compute the RHS term
figure
plot(sol-uu)
% full([sol uu])
% end

