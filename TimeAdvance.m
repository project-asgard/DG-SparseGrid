function f = TimeAdvance(A,f,t,dt,compression,Deg,pde,HASHInv)
%-------------------------------------------------
% Time Advance Method Input: Matrix:: A
%        Vector:: f Time Step:: dt
% Output: Vector:: f
%-------------------------------------------------

if pde.implicit   
    compression = 1;  
    %f = backward_euler(A,f,t,dt,compression,Deg,pde,HASHInv);
    f = crank_nicolson(A,f,t,dt,compression,Deg,pde,HASHInv);    
else   
    f = RungeKutta3(A,f,t,dt,compression,Deg,pde,HASHInv);   
end

end

function fval = RungeKutta3(A,f,t,dt,compression,Deg,pde,HASHInv)
%----------------------------------
% 3-rd Order Kutta Method
%----------------------------------

c2 = 1/2; c3 = 1;
source1 = source_vector(HASHInv,pde,t);
source2 = source_vector(HASHInv,pde,t+c2*dt);
source3 = source_vector(HASHInv,pde,t+c3*dt);

a21 = 1/2;
a31 = -1;
a32 = 2;
b1 = 1/6;
b2 = 2/3;
b3 = 1/6;

k_1 = ApplyA(A,f,compression,Deg,pde)   + source1;
y_1 = f + dt*a21*k_1;
k_2 = ApplyA(A,y_1,compression,Deg,pde) + source2;
y_2 = f+ dt*(a31*k_1+a32*k_2);
k_3 = ApplyA(A,y_2,compression,Deg,pde) + source3;

fval = f + dt*(b1*k_1+b2*k_2+b3*k_3);

end

function f1 = backward_euler(A,f0,t,dt,compression,Deg,pde,HASHInv)
%----------------------------------
% Backward Euler (First Order Implicit Time Advance)
%----------------------------------

s1 = source_vector(HASHInv,pde,t+dt);

[~,AMat] = ApplyA(A,f0,1,Deg,pde);

I = eye(numel(diag(AMat)));
AA = I - dt*AMat;

b = f0 + dt*s1;

f1 = AA\b; % Solve at each timestep

end

function f1 = crank_nicolson(A,f0,t,dt,compression,Deg,pde,HASHInv)
%----------------------------------
% Crank Nicolson (Second Order Implicit Time Advance)
%----------------------------------

s0 = source_vector(HASHInv,pde,t);
s1 = source_vector(HASHInv,pde,t+dt);

[~,AMat] = ApplyA(A,f0,1,Deg,pde);

I = eye(numel(diag(AMat)));
AA = 2*I - dt*AMat;

b = 2*f0 + dt*AMat*f0 + dt*(s0+s1);

f1 = AA\b; % Solve at each timestep

end

function [ftmp,A] = ApplyA(A_Data,f,compression,Deg,pde)

%-----------------------------------
% Multiply Matrix A by Vector f
%-----------------------------------
dof = size(f,1);
ftmp=sparse(dof,1);
Dim = pde.Dim;

if compression == 0
    
    % Explicitly construct the full A matrix (largest memory)
    
    dof = numel(A_Data.element_global_row_index);
    
    dofCnt = 1;
    conCnt = 1;
    
    A = zeros(dof,dof);
    
    for i=1:dof
        
        nConnected = A_Data.element_n_connected(i);
        
        for j=1:nConnected
            
            Index_I1 = A_Data.element_local_1_index(dofCnt);
            Index_I2 = A_Data.element_local_2_index(dofCnt);
            Index_J1 = A_Data.connected_local_1_index(conCnt);
            Index_J2 = A_Data.connected_local_2_index(conCnt);
            
            IndexI = A_Data.element_global_row_index(dofCnt);
            IndexJ = A_Data.connected_global_col_index(conCnt);
            
            % Apply term1 v.d_dx (vMassV . GradX)
            
            tmpA = A_Data.vMassV(Index_I1,Index_J1);
            tmpB = A_Data.GradX(Index_I2,Index_J2);
            
            A(IndexI,IndexJ) = A(IndexI,IndexJ) + tmpA*tmpB;
            
            % Apply term 2 E.d_dv (EMassX . GradV)
            
            tmpA = A_Data.GradV(Index_I1,Index_J1);
            tmpB = A_Data.EMassX(Index_I2,Index_J2);
            
            A(IndexI,IndexJ) = A(IndexI,IndexJ) + tmpA*tmpB;
            
            conCnt = conCnt+1;
            
        end
        
        dofCnt = dofCnt + 1;
        
    end
    
    % Do the matrix-vector multiply
    
    
    
    ftmp = sparse(A*f);
    
    
elseif compression == 1
    
    % Explicitly construct the sparse A matrix
    
    dof = numel(A_Data.element_global_row_index);
    
    dofCnt = 1;
    conCnt = 1;
    
    A = sparse(dof,dof);
    
    for i=1:dof
        
        nConnected = A_Data.element_n_connected(i);
        
        for j=1:nConnected
            
            Index_I1 = A_Data.element_local_1_index(dofCnt);
            Index_I2 = A_Data.element_local_2_index(dofCnt);
            Index_J1 = A_Data.connected_local_1_index(conCnt);
            Index_J2 = A_Data.connected_local_2_index(conCnt);
            
            IndexI = A_Data.element_global_row_index(dofCnt);
            IndexJ = A_Data.connected_global_col_index(conCnt);
            
            % Apply term1 v.d_dx (vMassV . GradX)
            
            tmpA = A_Data.vMassV(Index_I1,Index_J1);
            tmpB = A_Data.GradX(Index_I2,Index_J2);
            
            A(IndexI,IndexJ) = A(IndexI,IndexJ) + tmpA*tmpB;
            
            % Apply term 2 E.d_dv (EMassX . GradV)
            
            tmpA = A_Data.GradV(Index_I1,Index_J1);
            tmpB = A_Data.EMassX(Index_I2,Index_J2);
            
            A(IndexI,IndexJ) = A(IndexI,IndexJ) + tmpA*tmpB;
            
            conCnt = conCnt+1;
            
        end
        
        dofCnt = dofCnt + 1;
        
    end
    
    % Do the matrix-vector multiply
    
    ftmp = A*f;
    
    
elseif compression == 2
    
    % Elementwise matrix-vector multipliction (no tensor product encoding),
    % i.e., tmpA and tmpB are scalars
    
    dof = numel(A_Data.element_global_row_index);
    
    dofCnt = 1;
    conCnt = 1;
    
    for i=1:dof
        
        nConnected = A_Data.element_n_connected(i);
        
        for j=1:nConnected
            
            Index_I1 = A_Data.element_local_1_index(dofCnt);
            Index_I2 = A_Data.element_local_2_index(dofCnt);
            Index_J1 = A_Data.connected_local_1_index(conCnt);
            Index_J2 = A_Data.connected_local_2_index(conCnt);
            
            IndexI = A_Data.element_global_row_index(dofCnt);
            IndexJ = A_Data.connected_global_col_index(conCnt);
            
            % Apply term1 v.d_dx (vMassV . GradX)
            
            tmpA = A_Data.vMassV(Index_I1,Index_J1);
            tmpB = A_Data.GradX(Index_I2,Index_J2);
            
            ftmp(IndexI)=ftmp(IndexI)+tmpA*tmpB*f(IndexJ);
            
            % Apply term 2 E.d_dv (EMassX . GradV)
            
            tmpA = A_Data.GradV(Index_I1,Index_J1);
            tmpB = A_Data.EMassX(Index_I2,Index_J2);
            
            ftmp(IndexI)=ftmp(IndexI)+tmpA*tmpB*f(IndexJ);
            
            conCnt = conCnt+1;
            
        end
        
        dofCnt = dofCnt + 1;
        
    end
    
    
elseif compression == 3
    
    % Tensor product encoding over Deg (A_encode),
    % i.e., tmpA and tmpB are Deg x Deg matricies
    
    % Note: here A_Data == A_encode and follows the A_encode data
    % structure.
    
    for i=1:size(A_Data,2)
        
        tmpA=A_Data{i}.A1;
        tmpB=A_Data{i}.A2;
        IndexI=A_Data{i}.IndexI;
        IndexJ=A_Data{i}.IndexJ;
        
        ftmp(IndexI)=ftmp(IndexI)+kron_mult2(tmpA,tmpB,f(IndexJ));
        
    end
    
elseif compression == 4
    
    % Tensor product encoding over Deg (A_Data),
    % i.e., tmpA and tmpB are Deg x Deg matricies
    
    nWork = numel(A_Data.element_global_row_index);
    
    workCnt = 1;
    conCnt = 1;
    
    for workItem=1:nWork
        
        nConnected = A_Data.element_n_connected(workItem);
        
        
%         element_idx1D_1 = A_Data.element_local_1_index(workCnt);
%         element_idx1D_2 = A_Data.element_local_2_index(workCnt);
        element_idx1D = A_Data.element_local_index(workCnt,:);
        
        % Expand out the local and global indicies for this compressed item
        
%         Index_I1 = zeros(Deg,1);
%         Index_I2 = zeros(Deg,1);
        
        globalRow = zeros(Deg^2,1);
        degCnt1 = 1;
        degCnt2 = 1;
        for k1 = 1:Deg
%             Index_I1(degCnt1) = (element_idx1D_1-1)*Deg+k1;
%             Index_I2(degCnt1) = (element_idx1D_2-1)*Deg+k1;
            Index_I(degCnt1,:) = (element_idx1D-1)*Deg+k1;
            for k2 = 1:Deg
                globalRow(degCnt2) = Deg^2*(workItem-1)+Deg*(k1-1)+k2;
                degCnt2 = degCnt2 + 1;
            end
            degCnt1 = degCnt1 + 1;
        end
        
        for j=1:nConnected
            
%             connected_idx1D_1 = A_Data.connected_local_1_index(conCnt);
%             connected_idx1D_2 = A_Data.connected_local_2_index(conCnt);
            connected_idx1D = A_Data.connected_local_index(conCnt,:);
            connectedCol = A_Data.connected_global_col_index(conCnt);
            
            % Expand out the global col indicies for this compressed
            % connected item.
            
%             Index_J1 = zeros(Deg,1);
%             Index_J2 = zeros(Deg,1);
            globalCol = zeros(Deg^2,1);
            degCnt1 = 1;
            degCnt2 = 1;
            for kk1 = 1:Deg
%                 Index_J1(degCnt1) = (connected_idx1D_1-1)*Deg+kk1;
%                 Index_J2(degCnt1) = (connected_idx1D_2-1)*Deg+kk1;
                Index_J(degCnt1,:) = (connected_idx1D-1)*Deg+kk1;
                for kk2 = 1:Deg
                    globalCol(degCnt2) = Deg^2*(connectedCol-1)+Deg*(kk1-1)+kk2;
                    degCnt2 = degCnt2 + 1;
                end
                degCnt1 = degCnt1 + 1;
            end
            
            % Apply term1 v.d_dx (vMassV . GradX)
            if Dim == 2 
%             tmpA = A_Data.vMassV(Index_I1,Index_J1);
%             tmpB = A_Data.GradX(Index_I2,Index_J2);
            tmpA = A_Data.vMassV(Index_I(:,1),Index_J(:,1));
            tmpB = A_Data.GradX(Index_I(:,2),Index_J(:,2));
            end
            
            
            ftmp(globalRow)=ftmp(globalRow)+kron_mult2(tmpA,tmpB,f(globalCol));
            
            % Apply term 2 E.d_dv (EMassX . GradV)
            if Dim == 2 
%             tmpA = A_Data.GradV(Index_I1,Index_J1);
%             tmpB = A_Data.EMassX(Index_I2,Index_J2);
            tmpA = A_Data.GradV(Index_I(:,1),Index_J(:,1));
            tmpB = A_Data.EMassX(Index_I(:,2),Index_J(:,2));
            end
            
            ftmp(globalRow)=ftmp(globalRow)+kron_mult2(tmpA,tmpB,f(globalCol));
            
            conCnt = conCnt+1;
            
        end
        
        workCnt = workCnt + 1;
        
    end
end

end
