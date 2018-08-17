% This code solves Maxwells equations
% curl (curl u)-k^2 u = f
% with Dirichlet boundary conditions
%
% Ref: P. Houston, I. Perugia, A. Schneebeli, D. Schotzau,
% Interior penalty method for the indenitetime-harmonic Maxwell equations
% 2005

% Full Grids with DG scaling basis
% Full Grids with Multi-wavelet basis
% Sparse Grids

clear
close all
format short e
global A_encode

addpath(genpath(pwd))

Lev = 2;
Deg = 3;
Lmax = 1;
pde = Maxwell1;

TypeSolver = 'cg';
TypeMesh = 'sparsegrids';%'fullgrids';%
TypeAssembleMatrix = 'elementwise';%'';%

% DG scheme as following:
% (curl U,curl V)_T-k^2(U,V)_T
% -<{curl U},[Vxn]>_F-sigma*<{curl V},[Uxn]>_F+alpha*<[Uxn],[Vxn]>_F
% +alpha_1*<[curl Uxn],[curl Vxn]>_F
% = (f,V)
% parameter for DG penalty
% alpha h^{-1}<[Uxn],[Vxn]>_F
alpha = (2*Deg*(Deg+2))*100;%1000*(Deg)*(Deg+1);
alpha1 = 0;%1e-2;%1e-1;%1e-4;
% IPDG symmetric or non-symmetric
% sigma = 1 --> symmetric
% sigma =-1 --> non-symmetric
sigma = 1;


% if use iterative solver
if isequal(TypeSolver,'cg') == 1
    MaxIter = 1000;
    Tol = 1e-3;
elseif isequal(TypeSolver,'gmres') == 1
    MaxIter = 1000;
    Tol = 1e-3;
elseif isequal(TypeSolver,'multigrid')==1
    % set parameters for multigrid
end

% Coefficient Matrix
Mat1D = Matrix4CurlCurlOperator(Lev,Deg,Lmax);
Rhs1D = Rhs4CurlCurl(Deg,Lev,Lmax,pde);

% Meshing
% combs(Dim_i,:) means the lev for ith dimension
if isequal(TypeMesh,'fullgrids') == 1
    val = [0:Lev];
    n3 = repmat(val,1,length(val)*length(val));
    n3 = n3(:);
    n2 = repmat(val,length(val),length(val));
    n2 = n2(:);
    n1 = repmat(val,length(val)*length(val),1);
    n1 = n1(:);
    
    combs = [n1,n2,n3];
elseif isequal(TypeMesh,'sparsegrids') == 1
    Dim = 3;
    combs = perm_leq(Dim,Lev);
end


% index for the grids
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

Dofs = Deg^3*count;
Dofs3 = 3*Dofs;

% clear A_encode
A_encode={};

% generate the combination of lev
if isequal(TypeAssembleMatrix,'elementwise') == 1
    for l1 = 1:size(combs,1)
        
        row1 = Deg*(2^combs(l1,1)-2^max(0,combs(l1,1)-1))+1:Deg*2^combs(l1,1);
        row2 = Deg*(2^combs(l1,2)-2^max(0,combs(l1,2)-1))+1:Deg*2^combs(l1,2);
        row3 = Deg*(2^combs(l1,3)-2^max(0,combs(l1,3)-1))+1:Deg*2^combs(l1,3);
        
        IndexI = [Deg^3*(combs_index(l1,1)-1)+1:Deg^3*combs_index(l1,2)];
        
        for l2 = 1:size(combs,1)
            
            col1 = Deg*(2^combs(l2,1)-2^max(0,combs(l2,1)-1))+1:Deg*2^combs(l2,1);
            col2 = Deg*(2^combs(l2,2)-2^max(0,combs(l2,2)-1))+1:Deg*2^combs(l2,2);
            col3 = Deg*(2^combs(l2,3)-2^max(0,combs(l2,3)-1))+1:Deg*2^combs(l2,3);
            
            IndexJ = [Deg^3*(combs_index(l2,1)-1)+1:Deg^3*combs_index(l2,2)];
            
            A_loc = AssembleAencodeMaxwell3(IndexI,IndexJ,Mat1D,pde.w2,Dofs,...
                row1,row2,row3,col1,col2,col3,alpha,alpha1);
            
            A_encode = [A_encode,A_loc];
            
            
        end
        
        
    end
else
    
    IndexI = [1:Deg*2^Lev];
    %     IndexJ = [1:Deg*2^Lev];
    row1 = IndexI;col1 = IndexI;
    row2 = IndexI;col2 = IndexI;
    row3 = IndexI;col3 = IndexI;
    IndexI = [1:Dofs];
    IndexJ = [1:Dofs];
    A_encode = AssembleAencodeMaxwell3(IndexI,IndexJ,Mat1D,pde.w2,Dofs,...
        row1,row2,row3,col1,col2,col3,alpha,alpha1);
end

% Mat = Aencode2Matrix(A_encode,Dofs3);

% return
% Assemble the RHS term for Sparse Grids
ff = sparse(Dofs3,1);
if isequal(TypeAssembleMatrix,'elementwise') == 1
    for l1 = 1:size(combs,1)
        
        row1 = Deg*(2^combs(l1,1)-2^max(0,combs(l1,1)-1))+1:Deg*2^combs(l1,1);
        row2 = Deg*(2^combs(l1,2)-2^max(0,combs(l1,2)-1))+1:Deg*2^combs(l1,2);
        row3 = Deg*(2^combs(l1,3)-2^max(0,combs(l1,3)-1))+1:Deg*2^combs(l1,3);
        
        IndexI = [Deg^3*(combs_index(l1,1)-1)+1:Deg^3*combs_index(l1,2)];
        
        ff(IndexI) = kron(kron(Rhs1D.f1x(row1),Rhs1D.f1y(row2)),Rhs1D.f1z(row3));
        ff(Dofs+IndexI) = kron(kron(Rhs1D.f2x(row1),Rhs1D.f2y(row2)),Rhs1D.f2z(row3));
        ff(2*Dofs+IndexI) = kron(kron(Rhs1D.f3x(row1),Rhs1D.f3y(row2)),Rhs1D.f3z(row3));
    end
else
    ff =[kron(kron(Rhs1D.f1x,Rhs1D.f1y),Rhs1D.f1z);...
        kron(kron(Rhs1D.f2x,Rhs1D.f2y),Rhs1D.f2z);...
        kron(kron(Rhs1D.f3x,Rhs1D.f3y),Rhs1D.f3z)];
        
end

ff = ff*(2*pi^2-pde.w2);

% Tol = Tol*norm(ff);
x = zeros(Dofs3,1);
uu = ff/(2*pi^2-pde.w2);

sol = cg(x,ff,MaxIter,Tol);
[max(abs(sol-uu)) norm(sol-uu)]

[sol_CG,flag2,rr2,iter2,rv_CG] = pcg(@afun,ff,Tol,MaxIter);
semilogy(0:length(rv_CG)-1,rv_CG/norm(ff),'-o','linewidth',2,'markersize',8);
hold on;
[max(abs(sol_CG-uu)) norm(sol_CG-uu)]

[sol_GMRES,fl0,rr0,it0,rv_RMRES]= gmres(@afun,ff,10,Tol,MaxIter);
semilogy(0:length(rv_RMRES)-1,rv_RMRES/norm(ff),'-o','linewidth',2,'markersize',8);
legend({'CG','GMRES'},'fontsize',20)
xlabel('Iteration number');
ylabel('Relative residual');
[max(abs(sol_GMRES-uu)) norm(sol_GMRES-uu)]



