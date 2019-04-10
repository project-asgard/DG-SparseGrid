function f = TimeAdvance(pde,runTimeOpts,A_data,f,t,dt,deg,HASHInv,Vmax,Emax)
%-------------------------------------------------
% Time Advance Method Input: Matrix:: A
%        Vector:: f Time Step:: dt
% Output: Vector:: f
%-------------------------------------------------

if runTimeOpts.implicit
    %f = backward_euler(pde,runTimeOpts,A_data,f,t,dt,deg,HASHInv,Vmax,Emax);
    f = crank_nicolson(pde,runTimeOpts,A_data,f,t,dt,deg,HASHInv,Vmax,Emax);
else
    f = RungeKutta3(pde,runTimeOpts,A_data,f,t,dt,deg,HASHInv,Vmax,Emax);
end

end

function fval = RungeKutta3(pde,runTimeOpts,A_data,f,t,dt,deg,HASHInv,Vmax,Emax)
%----------------------------------
% 3-rd Order Kutta Method
%----------------------------------

%%
% Sources
c2 = 1/2; c3 = 1;
source1 = source_vector(HASHInv,pde,t);
source2 = source_vector(HASHInv,pde,t+c2*dt);
source3 = source_vector(HASHInv,pde,t+c3*dt);

%%
% Inhomogeneous dirichlet boundary conditions
bc1 = getBoundaryCondition1(pde,HASHInv,t); 
bc2 = getBoundaryCondition1(pde,HASHInv,t+c2*dt); 
bc3 = getBoundaryCondition1(pde,HASHInv,t+c3*dt); 

a21 = 1/2; a31 = -1; a32 = 2;
b1 = 1/6; b2 = 2/3; b3 = 1/6;

k_1 = ApplyA(pde,runTimeOpts,A_data,f,deg,Vmax,Emax)   + source1;% - bc1;
y_1 = f + dt*a21*k_1;
k_2 = ApplyA(pde,runTimeOpts,A_data,y_1,deg,Vmax,Emax) + source2;% - bc2;
y_2 = f+ dt*(a31*k_1+a32*k_2);
k_3 = ApplyA(pde,runTimeOpts,A_data,y_2,deg,Vmax,Emax) + source3;% - bc3;

fval = f + dt*(b1*k_1+b2*k_2+b3*k_3);

end

function f1 = backward_euler(pde,runTimeOpts,A_data,f0,t,dt,deg,HASHInv,Vmax,Emax)
%----------------------------------
% Backward Euler (First Order Implicit Time Advance)
%----------------------------------

s1 = source_vector(HASHInv,pde,t+dt);
bc1 = getBoundaryCondition1(pde,HASHInv,t+dt);

[~,AMat] = ApplyA(pde,runTimeOpts,A_data,f0,deg);

I = eye(numel(diag(AMat)));
AA = I - dt*AMat;

b = f0 + dt*s1 - dt*bc1;

f1 = AA\b; % Solve at each timestep

end

function f1 = crank_nicolson(pde,runTimeOpts,A_data,f0,t,dt,deg,HASHInv,Vmax,Emax)
%----------------------------------
% Crank Nicolson (Second Order Implicit Time Advance)
%----------------------------------

s0 = source_vector(HASHInv,pde,t);
s1 = source_vector(HASHInv,pde,t+dt);

bc0 = getBoundaryCondition1(pde,HASHInv,t);
bc1 = getBoundaryCondition1(pde,HASHInv,t+dt);

[~,AMat] = ApplyA(pde,runTimeOpts,A_data,f0,deg);

I = eye(numel(diag(AMat)));
AA = 2*I - dt*AMat;

b = 2*f0 + dt*AMat*f0 + dt*(s0+s1) - dt*(bc0+bc1);

f1 = AA\b; % Solve at each timestep

end

function [ftmp,A] = ApplyA(pde,runTimeOpts,A_data,f,deg,Vmax,Emax)

%-----------------------------------
% Multiply Matrix A by Vector f
%-----------------------------------
dof = size(f,1);
use_sparse_ftmp = 0;
if (use_sparse_ftmp),
  ftmp=sparse(dof,1);
else
  ftmp = zeros(dof,1);
end;
use_kronmultd = 1;

nTerms = numel(pde.terms);
nDims = numel(pde.dimensions);

dimensions = pde.dimensions;

useConnectivity = runTimeOpts.useConnectivity;

if runTimeOpts.compression == 0
    
    % Explicitly construct the full A matrix (largest memory)
    
    % Don't need this so removed. runTimeOpts.compression == 1 is fine.
    
    error('runTimeOpts.compression == 0 no longer valid, use runTimeOpts.compression == 4');
    
elseif runTimeOpts.compression == 1
    
    % Explicitly construct the sparse A matrix
    
    % Don't need since matrix construction for implicit is now done more
    % logically within the runTimeOpts.compression = 4.
    
    error('runTimeOpts.compression == 1 no longer valid, use runTimeOpts.compression == 4');
    
elseif runTimeOpts.compression == 2
    
    % Not used.
    
    error('runTimeOpts.compression == 2 no longer valid, use runTimeOpts.compression == 4');
    
elseif runTimeOpts.compression == 3
    
    % Tensor product encoding over Deg (A_encode),
    % i.e., tmpA and tmpB are Deg x Deg matricies
    
    % Note: here A_Data == A_encode and follows the A_encode data
    % structure.
    
    for i=1:size(A_data,2)
        
        tmpA=A_data{i}.A1;
        tmpB=A_data{i}.A2;
        IndexI=A_data{i}.IndexI;
        IndexJ=A_data{i}.IndexJ;
        
        if (use_kronmultd)
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
    
elseif runTimeOpts.compression == 4
    
    %%
    % Tensor product encoding over DOF within an element, i.e., over "deg" (A_Data),
    % i.e., tmpA and tmpB are deg_1 x deg_2 x deg_D matricies
    
    nWork = numel(A_data.element_global_row_index);
    
    conCnt = 1;
    
    ftmpA = ftmp;
    
    elementDOF = deg^nDims;
    
    implicit = runTimeOpts.implicit;

    if implicit
        totalDOF = nWork * elementDOF;
        A = sparse(totalDOF,totalDOF);
    end
      
    for workItem=1:nWork
        
        %disp([num2str(workItem) ' of ' num2str(nWork)]);
        
        if useConnectivity
            nConnected = A_data.element_n_connected(workItem);
        else
            nConnected = nWork; % Simply assume all are connected.
        end
        
        for d=1:nDims
            element_idx1D_D{d} = A_data.element_local_index_D{d}(workItem);
        end
        
        % Expand out the local and global indicies for this compressed item
        
        %%
        % TODO : add dimension dependent deg, something like ...
        % elementDOF = 1;
        % for d=1:nDims
        %     elementDOF = elementDOF * dimensions{d}.deg;
        % end
        
        globalRow = elementDOF*(workItem-1) + [1:elementDOF]';
        
        for d=1:nDims
            myDeg = dimensions{d}.deg;
            Index_I{d} = (element_idx1D_D{d}-1)*myDeg + [1:myDeg]';
        end
        
        for j=1:nConnected
            
            for d=1:nDims
                if useConnectivity
                    connected_idx1D_D{d} = A_data.connected_local_index_D{d}(conCnt);
                else
                    connected_idx1D_D{d} = A_data.element_local_index_D{d}(j);
                end
            end
            
            if useConnectivity
                connectedCol = A_data.connected_global_col_index(conCnt);
            else
                connectedCol = j;
            end
            
            % Expand out the global col indicies for this compressed
            % connected item.
            
            % NOTE : if we go to p-adaptivity then we will need 
            % a connected element DOF (connElementDOF) or the like.
            
            globalCol = elementDOF*(connectedCol-1) + [1:elementDOF]';
            
            for d=1:nDims
                myDeg = dimensions{d}.deg;
                Index_J{d} = (connected_idx1D_D{d}-1)*myDeg + [1:myDeg]';
            end
            
            %%
            % Apply operator matrices to present state using the pde spec
            % Y = A * X
            % where A is tensor product encoded.
            
            for t=1:nTerms
            
                %%
                % Construct the list of matrices for the kron_mult for this
                % operator (which has dimension many entries).
                clear kronMatList;
                for d=1:nDims                    
                    idx_i = Index_I{d}; 
                    idx_j = Index_J{d};
                    tmp = pde.terms{t}{d}.coeff_mat;
                    kronMatList{d} = tmp(idx_i,idx_j); % List of tmpA, tmpB, ... tmpD used in kron_mult
                end
                
                if implicit
                
                    %%
                    % Apply krond to return A (implicit time advance)
                    
                    A(globalRow,globalCol) = A(globalRow,globalCol) + krond(nDims,kronMatList);
                    
                else
                    
                    %%
                    % Apply kron_mult to return A*Y (explicit time advance)
                    X = f(globalCol);
                    if use_kronmultd
                        Y = kron_multd(nDims,kronMatList,X);
                    else
                        Y = kron_multd_full(nDims,kronMatList,X);
                    end

                    use_globalRow = 0;
                    if (use_globalRow),
                      ftmpA(globalRow) = ftmpA(globalRow) + Y;
                    else
                      % ------------------------------------------------------
                      % globalRow = elementDOF*(workItem-1) + [1:elementDOF]';
                      % ------------------------------------------------------
                      i1 = elementDOF*(workItem-1) + 1;
                      i2 = elementDOF*(workItem-1) + elementDOF;
                      ftmpA(i1:i2) = ftmpA(i1:i2) + Y;
                     end;
                    
                end
                
                
            end
            
            
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
                
    end
end

end
