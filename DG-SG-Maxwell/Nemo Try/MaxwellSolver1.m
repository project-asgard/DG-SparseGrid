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

%% Assembling the global Maxwell Matrix
for i = 1:Dof_Hash
    ll = InvHash{i};

    % Lev, Cell, Index from Hash
    iLev1 = ll(1);iLev2 = ll(2);iLev3 = ll(3);
    iCel1 = ll(4);iCel2 = ll(5);iCel3 = ll(6);
    iInd1 = ll(7);iInd2 = ll(8);iInd3 = ll(9);
    
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
            tmp_matB = -kron(tmp_matB12,tmp_matB3);
            
            
            II = Deg^3*(i-1)+1:Deg^3*i;
            JJ = Deg^3*(tmp-1)+1:Deg^3*tmp;
            [II,JJ] = meshgrid(JJ,II);
            
            Iu = [II,II+Dofs1,II+Dofs3,II+Dofs3+Dofs1];
            Iv = [JJ+Dofs3+Dofs1,JJ+Dofs3,JJ+Dofs1,JJ];
            Aij = [tmp_matB/(eps*mu),...
                  -tmp_matB/(eps*mu),...
                   tmp_matE,...
                  -tmp_matE];

            MaxMat = MaxMat+sparse(Iu,Iv,Aij,Dofs,Dofs);
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
            tmp_matE = -kron(kron(tmp_matE1,tmp_matE2),tmp_matE3);
            tmp_matB2 = Grad2X(Deg*(iInd2-1)+1:Deg*iInd2,Deg*(jIndex1D(j)-1)+1:Deg*jIndex1D(j));
            tmp_matB1 = eye(Deg);tmp_matB3=eye(Deg);
            tmp_matB = kron(kron(tmp_matB1,tmp_matB2),tmp_matB3);
            
            II = Deg^3*(i-1)+1:Deg^3*i;
            JJ = Deg^3*(tmp-1)+1:Deg^3*tmp;
            [II,JJ] = meshgrid(JJ,II);

            Iu = [II,II+2*Dofs1,II+Dofs3,II+Dofs3+2*Dofs1];
            Iv = [JJ+Dofs3+2*Dofs1,JJ+Dofs3,JJ+2*Dofs1,JJ];
            Aij =  [tmp_matB/(eps*mu),...
                   -tmp_matB/(eps*mu),...
                    tmp_matE,...
                   -tmp_matE];

            MaxMat = MaxMat+sparse(Iu,Iv,Aij,Dofs,Dofs);
            
            
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
            tmp_matB = -kron(tmp_matB1,tmp_matB23);
            
            II = Deg^3*(i-1)+1:Deg^3*i;
            JJ = Deg^3*(tmp-1)+1:Deg^3*tmp;
            [II,JJ] = meshgrid(JJ,II);

            Iu = [II+Dofs1,II+2*Dofs1,II+Dofs3+Dofs1,II+Dofs3+2*Dofs1];
            Iv = [JJ+Dofs3+2*Dofs1,JJ+Dofs3+Dofs1,JJ+2*Dofs1,JJ+Dofs1];
            Aij = [tmp_matB/(eps*mu),...
                  -tmp_matB/(eps*mu),...
                   tmp_matE,...
                  -tmp_matE];

            MaxMat = MaxMat+sparse(Iu,Iv,Aij,Dofs,Dofs);
            
        end
    end
    
    
end

%% Time advance for solving Maxwell equation
MaxRhs = -[Rhs/(eps);zeros(Dofs3,1)];
MaxSol = [E0;B0];

time=0;
for T=1:MaxT
    time=time+dt;

    bbb=MaxRhs*sin(omega*time);

    MaxSol = TimeAdvance(-MaxMat,MaxSol,dt,bbb);
    
end

Eh = MaxSol(1:Dofs3);
Bh = MaxSol(Dofs3+1:end);
