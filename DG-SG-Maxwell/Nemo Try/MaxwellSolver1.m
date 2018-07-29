function [Eh,Bh] = MaxwellSolver1(Lev,Deg,Hash,InvHash,Con1D,Grad1X,Grad2X,eps,mu,omega,dt,MaxT,...
                                  Rhs,E0,B0)
%=====================================================
% Maxwell Solver on [0,1]^3
% We use operator (du/dx,v) to construct (curl u,v)
% Note: The 3D matrix is not symmetric
%    [ 0,   B12, B13]
% Bs=[-B12,  0,  B23]
%    [-B13,-B23,  0 ]
%    [ 0,   E12, E23]
% Es=[-E12,  0,  E23]
%    [-E13,-E23,  0 ]
% MaxE=[-es/eps,0]
%        [  0   Es/(eps*mu)]
% MaxMat=[ Bs      0       ]
% Then time stepping for
% d [Eh]        [Eh] [MaxE]
%---    =MaxMat*    -
% dt[Bh]        [Bh] [0   ]
% Input:
%       Lev,Deg,Hash,InvHash,Con1D,eps,mu,omega
%       dt,MaxT: time stepping; steps
%       GradX: Grad Matrix
%       Rhs: bs
%       E0,B0: initials for E and B
% Output:
%       Eh
%       Bh
%=====================================================
% Assembling MaxMat 

Dof_Hash = size(InvHash,2);
Dofs1 = Dof_Hash*Deg^3;
Dofs3 = 3*Dofs1;
Dofs = 2*Dofs3;

MaxMat = sparse(Dofs,Dofs);
MaxRhs = sparse (Dofs,1);
MaxSol = sparse (Dofs,1);
ComMat1E = sparse(Dofs1,Dofs1);
ComMat2E = sparse(Dofs1,Dofs1);
ComMat3E = sparse(Dofs1,Dofs1);
ComMat1B = sparse(Dofs1,Dofs1);
ComMat2B = sparse(Dofs1,Dofs1);
ComMat3B = sparse(Dofs1,Dofs1);
%% Assembling the global Maxwell Matrix
for i = 1:Dof_Hash
    ll = InvHash{i};

    % Lev, Cell, Index from Hash
    iLev1 = ll(1);iLev2 = ll(2);iLev3 = ll(3);
    iCel1 = ll(4);iCel2 = ll(5);iCel3 = ll(6);
    iInd1 = ll(7);iInd2 = ll(8);iInd3 = ll(9);
    
    Mat1E=sparse(Dofs1,Dofs1);
    Mat2E=sparse(Dofs1,Dofs1);
    Mat3E=sparse(Dofs1,Dofs1);
    Mat1B=sparse(Dofs1,Dofs1);
    Mat2B=sparse(Dofs1,Dofs1);
    Mat3B=sparse(Dofs1,Dofs1);
    %****************************
    % consider A12 and A21=A12'
    % I*I*T
    %****************************
    jIndex1D = Con1D{iInd3};
    
    for j = 1:size(jIndex1D,2)
        
        jLev1=iLev1;jLev2=iLev2;jLev3=ceil(log2(jIndex1D(j)));
        jCel1=iCel1;jCel2=iCel2;jCel3=max(jIndex1D(j)-1-2^(jLev3-1),0);
        
        
        if jLev1+jLev2+jLev3<=Lev
            jkey = [jLev1,jLev2,jLev3,jCel1,jCel2,jCel3];

            tmp = Hash.(sprintf('i%g_',jkey));
            tmp_matE3 = Grad1X(Deg*(iInd3-1)+1:Deg*iInd3,Deg*(jIndex1D(j)-1)+1:Deg*jIndex1D(j));
            tmp_matE12 = speye(Deg^2);
            tmp_matE = kron(tmp_matE12,tmp_matE3);
            tmp_matB3 = Grad2X(Deg*(iInd3-1)+1:Deg*iInd3,Deg*(jIndex1D(j)-1)+1:Deg*jIndex1D(j));
            tmp_matB12 = speye(Deg^2);
            tmp_matB = kron(tmp_matB12,tmp_matB3);
            
            
            II = Deg^3*(i-1)+1:Deg^3*i;
            JJ = Deg^3*(tmp-1)+1:Deg^3*tmp;
            [II,JJ] = meshgrid(JJ,II);
            
            Mat1B=Mat1B+sparse(II,JJ,tmp_matB,Dofs1,Dofs1);
            Mat1E=Mat1E+sparse(II,JJ,tmp_matE,Dofs1,Dofs1);
        end
        
        
    end

    %****************************
    % consider A13 and A31=A13'
    % I*T*I
    %****************************
    jIndex1D = Con1D{iInd2};

    for j = 1:size(jIndex1D,2)
        
        jLev1=iLev1;jLev2=ceil(log2(jIndex1D(j)));jLev3=iLev3;
        jCel1=iCel1;jCel2=max(jIndex1D(j)-1-2^(jLev2-1),0);jCel3=iCel3;
        
        if jLev1+jLev2+jLev3<=Lev
            jkey = [jLev1,jLev2,jLev3,jCel1,jCel2,jCel3];
            % tmp is not contigent
            tmp = Hash.(sprintf('i%g_',jkey));
            tmp_matE2 = Grad1X(Deg*(iInd2-1)+1:Deg*iInd2,Deg*(jIndex1D(j)-1)+1:Deg*jIndex1D(j));
            tmp_matE1 = eye(Deg);tmp_matE3=eye(Deg);
            tmp_matE = kron(kron(tmp_matE1,tmp_matE2),tmp_matE3);
            tmp_matB2 = Grad2X(Deg*(iInd2-1)+1:Deg*iInd2,Deg*(jIndex1D(j)-1)+1:Deg*jIndex1D(j));
            tmp_matB1 = eye(Deg);tmp_matB3=eye(Deg);
            tmp_matB = kron(kron(tmp_matB1,tmp_matB2),tmp_matB3);
            
            II = Deg^3*(i-1)+1:Deg^3*i;
            JJ = Deg^3*(tmp-1)+1:Deg^3*tmp;
            [II,JJ] = meshgrid(JJ,II);

            Mat2B=Mat2B+sparse(II,JJ,tmp_matB,Dofs1,Dofs1);
            Mat2E=Mat2E+sparse(II,JJ,tmp_matE,Dofs1,Dofs1);
            
            
        end
    end
    
    %****************************
    % consider A23 and A32=A23'
    % T*I*I
    %****************************
    jIndex1D = Con1D{iInd1};

    for j = 1:size(jIndex1D,2)
        
        jLev1=ceil(log2(jIndex1D(j)));jLev2=iLev2;jLev3=iLev3;
        jCel1=max(jIndex1D(j)-1-2^(jLev1-1),0);jCel2=iCel2;jCel3=iCel3;
        
        if jLev1+jLev2+jLev3<=Lev
            jkey = [jLev1,jLev2,jLev3,jCel1,jCel2,jCel3];
            % tmp is not contigent
            tmp = Hash.(sprintf('i%g_',jkey));
            tmp_matE1 = Grad1X(Deg*(iInd1-1)+1:Deg*iInd1,Deg*(jIndex1D(j)-1)+1:Deg*jIndex1D(j));
            tmp_matE23 = speye(Deg^2);
            tmp_matE = kron(tmp_matE1,tmp_matE23);
            tmp_matB1 = Grad2X(Deg*(iInd1-1)+1:Deg*iInd1,Deg*(jIndex1D(j)-1)+1:Deg*jIndex1D(j));
            tmp_matB23 = speye(Deg^2);
            tmp_matB = kron(tmp_matB1,tmp_matB23);
            
            II = Deg^3*(i-1)+1:Deg^3*i;
            JJ = Deg^3*(tmp-1)+1:Deg^3*tmp;
            [II,JJ] = meshgrid(JJ,II);

            Mat3B=Mat3B+sparse(II,JJ,tmp_matB,Dofs1,Dofs1);
            Mat3E=Mat3E+sparse(II,JJ,tmp_matE,Dofs1,Dofs1);
            
        end
    end
    
        ComMat1E=ComMat1E+Mat1E;
        ComMat2E=ComMat2E+Mat2E;
        ComMat3E=ComMat3E+Mat3E;
        ComMat1B=ComMat1B+Mat1B;
        ComMat2B=ComMat2B+Mat2B;
        ComMat3B=ComMat3B+Mat3B;
        
end
zero=sparse(Dofs1,Dofs1);

MaxMat=[zero,zero,zero,zero,ComMat1B/(eps*mu),-ComMat2B/(eps*mu);
        zero,zero,zero,-ComMat1B/(eps*mu),zero,ComMat3B/(eps*mu);
        zero,zero,zero,ComMat2B/(eps*mu),-ComMat3B/(eps*mu),zero;
        zero,-ComMat1E,ComMat2E,zero,zero,zero;
        ComMat1E,zero,-ComMat3E,zero,zero,zero;
        -ComMat2E,ComMat3E,zero,zero,zero,zero;];
opts.tol=1e-4;
e=eigs(MaxMat,1,'lm',opts)
s=norm(MaxMat,'fro')
%% Time advance for solving Maxwell equation
MaxRhs = -[Rhs/(eps);zeros(Dofs3,1)];
MaxSol = [E0;B0];

time=0;
for T=1:MaxT
    time=time+dt;

    bbb=MaxRhs*sin(omega*time);

    MaxSol = TimeAdvance(MaxMat,MaxSol,dt,bbb);
    
end

Eh = MaxSol(1:Dofs3);
Bh = MaxSol(Dofs3+1:end);
