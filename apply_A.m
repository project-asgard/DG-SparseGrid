function [ftmp, A, ALHS] = apply_A (pde, opts, A_data, f, deg, Vmax, Emax)

%-----------------------------------
% Multiply Matrix A by Vector f
%-----------------------------------
dof = size(f,1);
use_sparse_ftmp = 0;
if (use_sparse_ftmp)
    ftmp=sparse(dof,1);
else
    ftmp = zeros(dof,1);
end
use_kronmultd = 1;

num_terms     = numel(pde.terms);
num_terms_LHS = numel(pde.termsLHS);
num_dims      = numel(pde.dimensions);

%%
% Tensor product encoding over DOF within an element, i.e., over "deg" (A_Data),
% i.e., tmpA and tmpB are deg_1 x deg_2 x deg_D matricies

nWork = numel(A_data.element_global_row_index);

conCnt = 1;

ftmpA = ftmp;

elementDOF = deg^num_dims;

implicit = opts.implicit;

totalDOF = nWork * elementDOF;
if opts.implicit
    A = zeros(totalDOF,totalDOF); % Only filled if implicit
end
ALHS = 0;
if num_terms_LHS > 0
    ALHS = zeros(totalDOF,totalDOF); % Only filled if non-identity LHS mass matrix
end

for workItem=1:nWork
          
    nConnected = nWork; % Simply assume all are connected. 
    
    for d=1:num_dims
        element_idx1D_D{d} = A_data.element_local_index_D{d}(workItem);
    end
    
    % Expand out the local and global indicies for this compressed item
    
    globalRow = elementDOF*(workItem-1) + [1:elementDOF]';
    
    for d=1:num_dims
        myDeg = pde.deg;
        Index_I{d} = (element_idx1D_D{d}-1)*myDeg + [1:myDeg]';
    end
    
    for j=1:nConnected
        
        for d=1:num_dims
            connected_idx1D_D{d} = A_data.element_local_index_D{d}(j);
        end
        
        connectedCol = j;

        % Expand out the global col indicies for this compressed
        % connected item.
        
        globalCol = elementDOF*(connectedCol-1) + [1:elementDOF]';
        
        for d=1:num_dims
            myDeg = pde.deg;
            Index_J{d} = (connected_idx1D_D{d}-1)*myDeg + [1:myDeg]';
        end
        
        %%
        % Apply operator matrices to present state using the pde spec
        % Y = A * X
        % where A is tensor product encoded.
        
        for t=1:num_terms
            
            %%
            % Construct the list of matrices for the kron_mult for this
            % operator (which has dimension many entries).
            for d=1:num_dims
                idx_i = Index_I{d};
                idx_j = Index_J{d};
                tmp = pde.terms{t}.terms_1D{d}.mat;
                kronMatList{d} = tmp(idx_i,idx_j); % List of tmpA, tmpB, ... tmpD used in kron_mult
            end
            
            if implicit
                
                %%
                % Apply krond to return A (implicit time advance)
                
                A(globalRow,globalCol) = A(globalRow,globalCol) + krond(num_dims,kronMatList);
                
            else
                
                %%
                % Apply kron_mult to return A*Y (explicit time advance)
                X = f(globalCol);
                if use_kronmultd
                    Y = kron_multd(num_dims,kronMatList,X);
                else
                    Y = kron_multd_full(num_dims,kronMatList,X);
                end
                
                use_globalRow = 0;
                if (use_globalRow)
                    ftmpA(globalRow) = ftmpA(globalRow) + Y;
                else
                    % ------------------------------------------------------
                    % globalRow = elementDOF*(workItem-1) + [1:elementDOF]';
                    % ------------------------------------------------------
                    i1 = elementDOF*(workItem-1) + 1;
                    i2 = elementDOF*(workItem-1) + elementDOF;
                    ftmpA(i1:i2) = ftmpA(i1:i2) + Y;
                end
                
            end
            
        end
        
        %%
        % Construct the mat list for a non-identity LHS mass matrix
        for t=1:num_terms_LHS
            for d=1:num_dims
                idx_i = Index_I{d};
                idx_j = Index_J{d};
                tmp = pde.termsLHS{t}.terms_1D{d}.mat;
                kronMatListLHS{d} = tmp(idx_i,idx_j); % List of tmpA, tmpB, ... tmpD used in kron_mult
            end
            
            %%
            % Apply krond to return A (recall this term requires inversion)
            
            ALHS(globalRow,globalCol) = ALHS(globalRow,globalCol) + krond(num_dims,kronMatListLHS);
            
        end
        
        
        %%
        % Overwrite previous approach with PDE spec approch
        ftmp = ftmpA;
        
        conCnt = conCnt+1;
        
    end
    
    assert(workItem==workItem);
    
end

end