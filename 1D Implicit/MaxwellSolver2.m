function [Bh,E1h,E2h] = MaxwellSolver2(Lev,Deg,Hash,InvHash,Con1D,GradX,...
                                 omega,dt,MaxT,...
                                 F_1D,E1,E2,B)%Gauss-Legendre method
%Still have some problem                        
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
% MaxRhs = -[zeros(Dofs,1);F_1D.b1;F_1D.b2];
MaxRhs=-[F_1D.b3;F_1D.b1;F_1D.b2];
MaxMat=sparse(Dofs*3,Dofs*3);
MaxSol = [B;E1;E2];
s1=inv(eye(Dofs*3)-sqrt(3)/6*MaxMat*dt);
s2=inv(eye(Dofs*3)+sqrt(3)/6*MaxMat*dt);
h1=inv(eye(Dofs*3)-1/4*MaxMat*dt-MaxMat*s1*(1/4-sqrt(3)/6)*(eye(Dofs*3)+sqrt(3)/6*MaxMat*dt)*dt);
h2=inv(eye(Dofs*3)-1/4*MaxMat*dt-MaxMat*s2*(1/4+sqrt(3)/6)*(eye(Dofs*3)-sqrt(3)/6*MaxMat*dt)*dt);
time=0;
%e=eigs(S,6,'largestabs','Tolerance',1e-2)
for T=1:MaxT
    time=time+dt;

%     b1=MaxRhs*sin(omega*(time+(1/2-1/6*sqrt(3))*dt));
%     b2=MaxRhs*sin(omega*(time+(1/2+1/6*sqrt(3))*dt));
%      b1=MaxRhs*(time-dt+(1/2-1/6*sqrt(3))*dt)^3;
%      b2=MaxRhs*(time-dt+(1/2+1/6*sqrt(3))*dt)^3;
     b1=MaxRhs*(time-dt+(1/2-1/6*sqrt(3))*dt);
     b2=MaxRhs*(time-dt+(1/2+1/6*sqrt(3))*dt);

    MaxSol = ImplicitTime2(s1,s2,h1,h2,MaxMat,MaxSol,dt,b1,b2);
    a=1;
    
    
    
end

Bh = MaxSol(1:Dofs);
E1h = MaxSol(Dofs+1:Dofs*2);
E2h = MaxSol(Dofs*2+1:end);
