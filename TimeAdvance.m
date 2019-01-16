function f = TimeAdvance(A,f,t,dt,compression,deg,pde,HASHInv,Vmax,Emax)
%-------------------------------------------------
% Time Advance Method Input: Matrix:: A
%        Vector:: f Time Step:: dt
% Output: Vector:: f
%-------------------------------------------------

if pde.implicit
    f = backward_euler(pde,A,f,t,dt,compression,deg,HASHInv,Vmax,Emax);
    %f = crank_nicolson(A,f,t,dt,compression,Deg,pde,HASHInv);
else
    f = RungeKutta3(pde,A,f,t,dt,compression,deg,HASHInv,Vmax,Emax);
end

end

function fval = RungeKutta3(pde,A,f,t,dt,compression,deg,HASHInv,Vmax,Emax)
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

k_1 = ApplyA(pde,A,f,compression,deg,Vmax,Emax)   + source1;
y_1 = f + dt*a21*k_1;
k_2 = ApplyA(pde,A,y_1,compression,deg,Vmax,Emax) + source2;
y_2 = f+ dt*(a31*k_1+a32*k_2);
k_3 = ApplyA(pde,A,y_2,compression,deg,Vmax,Emax) + source3;

fval = f + dt*(b1*k_1+b2*k_2+b3*k_3);

end

function f1 = backward_euler(pde,A,f0,t,dt,compression,deg,HASHInv,Vmax,Emax)
%----------------------------------
% Backward Euler (First Order Implicit Time Advance)
%----------------------------------

s1 = source_vector(HASHInv,pde,t+dt);

[~,AMat] = ApplyA(pde,A,f0,1,deg);

I = eye(numel(diag(AMat)));
AA = I - dt*AMat;

b = f0 + dt*s1;

f1 = AA\b; % Solve at each timestep

end

function f1 = crank_nicolson(pde,A,f0,t,dt,compression,deg,HASHInv,Vmax,Emax)
%----------------------------------
% Crank Nicolson (Second Order Implicit Time Advance)
%----------------------------------

s0 = source_vector(HASHInv,pde,t);
s1 = source_vector(HASHInv,pde,t+dt);

[~,AMat] = ApplyA(pde,A,f0,1,deg);

I = eye(numel(diag(AMat)));
AA = 2*I - dt*AMat;

b = 2*f0 + dt*AMat*f0 + dt*(s0+s1);

f1 = AA\b; % Solve at each timestep

end

function [ftmp,A] = ApplyA(pde,A_Data,f,compression,deg,Vmax,Emax)

%-----------------------------------
% Multiply Matrix A by Vector f
%-----------------------------------
dof = size(f,1);
ftmp=sparse(dof,1);
use_kronmult2 = 1;

nTerms = numel(pde.terms);
nDims = numel(pde.dimensions);

dimensions = pde.dimensions;

Identity = speye(dof,dof);

if compression == 0
    
    % Explicitly construct the full A matrix (largest memory)
    
    % Don't need this so removed. Compression == 1 is fine. 
    
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
            
            %%
            % TODO : this needs to be generalized to dim.
            
            Index_I{1} = Index_I1;
            Index_I{2} = Index_I2;
            
            Index_J{1} = Index_J1;
            Index_J{2} = Index_J2;
            
            %%
            % Apply operator matrices to present state using the pde spec
            % Y = A * X
            
            for t=1:nTerms
            
                %%
                % Construct the list of matrices for the kron_mult for this
                % operator (which has dimension many entries).
                Atmp = 1;
                for d=1:nDims                    
                    idx_i = Index_I{d}; 
                    idx_j = Index_J{d};
                    tmp = pde.terms{t}{d}.coeff_mat;
                    Atmp = Atmp * tmp(idx_i,idx_j);
                end
                
                %%
                % Construct global A matrix
                
                A(IndexI,IndexJ) = A(IndexI,IndexJ) + Atmp;

            end
            
%             %% 
%             % Apply each of the 2 terms in the equation 
%             
%             % Apply term1 v.d_dx (vMassV . GradX)
%             
%             tmpA = A_Data.vMassV(Index_I1,Index_J1);
%             tmpB = A_Data.GradX(Index_I2,Index_J2);
%             
%             A(IndexI,IndexJ) = A(IndexI,IndexJ) + tmpA*tmpB;
%             
%             % Apply term 2 E.d_dv (EMassX . GradV)
%             
%             tmpA = A_Data.GradV(Index_I1,Index_J1);
%             tmpB = A_Data.EMassX(Index_I2,Index_J2);
%             
%             A(IndexI,IndexJ) = A(IndexI,IndexJ) + tmpA*tmpB;
            
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
        
        if (use_kronmult2)
            ftmp(IndexI)=ftmp(IndexI)+kronmult2(tmpA,tmpB,f(IndexJ));
        else
            % [nrA,ncA] = size(tmpA);
            % [nrB,ncB] = size(tmpB);
            
            nrA = size(tmpA,1);
            ncA = size(tmpA,2);
            nrB = size(tmpB,1);
            ncB = size(tmpB,2);
            
            ftmp(IndexI)=ftmp(IndexI) + ...
                reshape(tmpB * reshape(f(IndexJ),ncB,ncA)*transpose(tmpA), nrB*nrA,1);
        end
        
    end
    
elseif compression == 4
    
    %%
    % Tensor product encoding over DOF within an element, i.e., over "deg" (A_Data),
    % i.e., tmpA and tmpB are deg_1 x deg_2 x deg_D matricies
    
    nWork = numel(A_Data.element_global_row_index);
    
%     workItem = 1;
    conCnt = 1;
    
    ftmpA = ftmp;
   
    for workItem=1:nWork
        
        nConnected = A_Data.element_n_connected(workItem);
        
        element_idx1D_1 = A_Data.element_local_1_index(workItem);
        element_idx1D_2 = A_Data.element_local_2_index(workItem);
        for d=1:nDims
            element_idx1D_D{d} = A_Data.element_local_index_D{d}(workItem);
        end
        
        % Expand out the local and global indicies for this compressed item
        
        Index_I1 = zeros(deg,1);
        Index_I2 = zeros(deg,1);
%       globalRow = zeros(Deg^2,1);
        globalRow = zeros(deg^nDims,1);
        degCnt1 = 1;
        degCnt2 = 1;
        for k1 = 1:deg
            Index_I1(degCnt1) = (element_idx1D_1-1)*deg+k1;
            Index_I2(degCnt1) = (element_idx1D_2-1)*deg+k1;
            for k2 = 1:deg
%               globalRow(degCnt2) = Deg^2*(workItem-1)+Deg*(k1-1)+k2;
                globalRow(degCnt2) = deg^nDims*(workItem-1)+deg*(k1-1)+k2;
                degCnt2 = degCnt2 + 1;
            end
            degCnt1 = degCnt1 + 1;
        end
        
        elementDOF = deg^nDims;
        
        %%
        % TODO : add dimension dependent deg, something like ...
        % elementDOF = 1;
        % for d=1:nDims
        %     elementDOF = elementDOF * dimensions{d}.deg;
        % end
        
        globalRowA = elementDOF*(workItem-1) + [1:elementDOF]';
        if nDims==2
            assert(norm(globalRowA - globalRow)==0);
        end
        globalRow = globalRowA;
        
        for d=1:nDims
            myDeg = dimensions{d}.deg;
            Index_I{d} = (element_idx1D_D{d}-1)*myDeg + [1:myDeg]';
        end
        if nDims==2
            assert(norm(Index_I{1}-Index_I1)==0);
            assert(norm(Index_I{2}-Index_I2)==0);
        end
        
        for j=1:nConnected
            
            connected_idx1D_1 = A_Data.connected_local_1_index(conCnt);
            connected_idx1D_2 = A_Data.connected_local_2_index(conCnt);
            for d=1:nDims
                connected_idx1D_D{d} = A_Data.connected_local_index_D{d}(conCnt);
            end
            
            connectedCol = A_Data.connected_global_col_index(conCnt);
            
            % Expand out the global col indicies for this compressed
            % connected item.
            
            Index_J1 = zeros(deg,1);
            Index_J2 = zeros(deg,1);
            globalCol = zeros(deg^2,1);
            degCnt1 = 1;
            degCnt2 = 1;
            for kk1 = 1:deg
                Index_J1(degCnt1) = (connected_idx1D_1-1)*deg+kk1;
                Index_J2(degCnt1) = (connected_idx1D_2-1)*deg+kk1;
                for kk2 = 1:deg
                    globalCol(degCnt2) = deg^2*(connectedCol-1)+deg*(kk1-1)+kk2;
                    degCnt2 = degCnt2 + 1;
                end
                degCnt1 = degCnt1 + 1;
            end
            
            % NOTE : if we go to p-adaptivity then we will need 
            % a connected element DOF (connElementDOF) or the like.
            
            globalColA = elementDOF*(connectedCol-1) + [1:elementDOF]';
            if nDims==2
                assert(norm(globalColA-globalCol)==0);
            end
            globalCol = globalColA;
            
            for d=1:nDims
                myDeg = dimensions{d}.deg;
                Index_J{d} = (connected_idx1D_D{d}-1)*myDeg + [1:myDeg]';
            end
            if nDims==2
                assert(norm(Index_J{1}-Index_J1)==0);
                assert(norm(Index_J{2}-Index_J2)==0);
            end
            
%             %%
%             % TODO : this needs to be generalized to dim.
%             
%             Index_I{1} = Index_I1;
%             Index_I{2} = Index_I2;
%             
%             Index_J{1} = Index_J1;
%             Index_J{2} = Index_J2;
            
            %%
            % Apply operator matrices to present state using the pde spec
            % Y = A * X
            % where A is tensor product encoded.
            
            for t=1:nTerms
            
                %%
                % Construct the list of matrices for the kron_mult for this
                % operator (which has dimension many entries).
                clear A;
                for d=1:nDims                    
                    idx_i = Index_I{d}; 
                    idx_j = Index_J{d};
                    tmp = pde.terms{t}{d}.coeff_mat;
                    A{d} = tmp(idx_i,idx_j); % List of tmpA, tmpB, ... tmpD used in kron_mult
                end
                
                %%
                % Apply kron_mult
                X = f(globalCol);
                if use_kronmult2
                    Y = kron_multd(nDims,A,X);
                else
                    Y = kron_multd_full(nDims,A,X);
                end
                ftmpA(globalRow) = ftmpA(globalRow) + Y;
            end
            
%             %%
%             % Apply operator matrices using original Vlasov hardwired
%             % approach
%             
%             % Apply term1 v.d_dx (vMassV . GradX)
%             
%             tmpA = A_Data.vMassV(Index_I1,Index_J1);
%             tmpB = A_Data.GradX(Index_I2,Index_J2);            
%             
%             if use_kronmult2
%                 ftmp(globalRow)=ftmp(globalRow)+kronmult2(tmpA,tmpB,f(globalCol));
%             else
%                 ftmp(globalRow)=ftmp(globalRow)+ ...
% 		      reshape(tmpB * reshape( f(globalCol),Deg,Deg)* transpose(tmpA),Deg*Deg,1);              
%             end          
%                          
%             % Apply term 2 E.d_dv (EMassX . GradV)
%             
%             tmpA = A_Data.GradV(Index_I1,Index_J1);
%             tmpB = A_Data.EMassX(Index_I2,Index_J2);
%             
%             if use_kronmult2
%                 ftmp(globalRow)=ftmp(globalRow)+kronmult2(tmpA,tmpB,f(globalCol));              
%             else
%                 ftmp(globalRow)=ftmp(globalRow)+ ...
% 		      reshape(tmpB * reshape( f(globalCol),Deg,Deg)* transpose(tmpA),Deg*Deg,1);               
%             end
            
            %%
            % Overwrite previous approach with PDE spec approch
            ftmp = ftmpA; 

                        
            % if CF, comment Line 360-386, if LF, uncomment them
%            % Apply term3 vMax*Mass x FluxX
%             
%             tmpA = Identity(Index_I1,Index_J1);
%             tmpB = Vmax*A_Data.FluxX(Index_I2,Index_J2);
%             
%             
%             if use_kronmult2
%                 ftmp(globalRow)=ftmp(globalRow)+kronmult2(tmpA,tmpB,f(globalCol));
%             else
%                 ftmp(globalRow)=ftmp(globalRow)+ ...
% 		      reshape(tmpB * reshape( f(globalCol),Deg,Deg)* transpose(tmpA),Deg*Deg,1);
%                 
%             end 
%             
%           % Apply term4 Emax*FluxV x Mass
%             
%             tmpA = -Emax*A_Data.FluxV(Index_I1,Index_J1);
%             tmpB = Identity(Index_I2,Index_J2);
%             
%             
%             if use_kronmult2
%                 ftmp(globalRow)=ftmp(globalRow)+kronmult2(tmpA,tmpB,f(globalCol));
%             else
%                 ftmp(globalRow)=ftmp(globalRow)+ ...
% 		      reshape(tmpB * reshape( f(globalCol),Deg,Deg)* transpose(tmpA),Deg*Deg,1);
%                 
%             end 
            
            conCnt = conCnt+1;
            
        end
        
        assert(workItem==workItem);
        
%         workItem = workItem + 1;
        
    end
end

end
