function [Bh,E1h,E2h] = MaxwellSolver1(Lev,Deg,Hash,InvHash,Con1D,GradX,...
                                 omega,dt,MaxT,...
                                 F_1D,E1,E2,B)%Trapezoidal rule
%=====================================================
% Maxwell Solver on [0,1]^3
% We use operator (du/dx,v) to construct (curl u,v)
% Note: The 3D matrix is not symmetric
%    [ 0,   A12, A13]
% As=[-A12,  0,  A23]
%    [-A13,-A23,  0 ]
%        [  0   As/(eps*nu)]
% MaxMat=[-As      0       ]
%
% MaxB=[-bs/eps,0]
% Then time stepping for
% d [Eh]        [Eh] [MaxB]
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
Dofs = Dof_Hash*Deg;

Mat = sparse(Dofs,Dofs);
MaxMat = sparse (Dofs*3,Dofs*3);
MaxRhs = sparse (Dofs*3,1);
MaxSol = sparse (Dofs*3,1);


%% Assembling the global Maxwell Matrix
for i = 1:Dof_Hash
    ll = InvHash{i};

    % Lev, Cell, Index from Hash
    iLev1 = ll(1);
    iCel1 = ll(2);
    iInd1 = ll(3);
  
    
    jIndex1D = Con1D{iInd1};

    for j = 1:size(jIndex1D,2)
        
        jLev1=ceil(log2(jIndex1D(j)));
        jCel1=max(jIndex1D(j)-1-2^(jLev1-1),0);
        
        if jLev1<=Lev
            jkey = [jLev1,jCel1];
            % tmp is not contigent
            tmp = Hash.(sprintf('i%g_',jkey));
            tmp_mat1 = GradX(Deg*(iInd1-1)+1:Deg*iInd1,Deg*(jIndex1D(j)-1)+1:Deg*jIndex1D(j));
            
            II = Deg*(i-1)+1:Deg*i;
            JJ = Deg*(tmp-1)+1:Deg*tmp;
            [II,JJ] = meshgrid(JJ,II);

            Mat = Mat+sparse(II,JJ,tmp_mat1,Dofs,Dofs);
            
        end
    end
    
    
end
zero=sparse(Dofs,Dofs);
MaxMat=-[zero,Mat,zero;
        Mat,zero,zero;
        zero,zero,zero];


%% Time advance for solving Maxwell equation
MaxRhs = -[zeros(Dofs,1);F_1D.b1;F_1D.b2];
MaxSol = [B;E1;E2];
s=eye(Dofs*3)-0.5*MaxMat*dt;
H=inv(s);
% % e=eigs(H*(eye(Dofs*3)+0.5*MaxMat*dt))
% % e=abs(e)
time=0;
for T=1:MaxT
    time=time+dt;

    bbb=MaxRhs*sin(omega*time);

    MaxSol = ImplicitTime1(H,MaxMat,MaxSol,dt,bbb);
    
end

Bh = MaxSol(1:Dofs);
E1h = MaxSol(Dofs+1:Dofs*2);
E2h = MaxSol(Dofs*2+1:end);




