function [ftmp, A, ALHS] = apply_A (pde,opts,A_data,f,deg,Vmax,Emax)

if opts.use_connectivity
    connectivity = pde.connectivity;
end

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

num_elem = numel(A_data.element_global_row_index);

ftmpA = ftmp;

elementDOF = deg^num_dims;

totalDOF = num_elem * elementDOF;
ALHS = 0;
A = 0;
if opts.build_A
    A = zeros(totalDOF,totalDOF); % Only filled if using hand coded implicit methods
end
if num_terms_LHS > 0
    ALHS = zeros(totalDOF,totalDOF); % Only filled if non-identity LHS mass matrix
end

if opts.build_A && ~opts.quiet; disp('Building A ...'); end
for elem=1:num_elem

    if opts.use_connectivity
        num_connected = numel(connectivity{elem});
    else
        num_connected = num_elem; % Simply assume all are connected.
    end
    
    for d=1:num_dims
        element_idx1D_D{d} = A_data.element_local_index_D{d}(elem);
    end
    
    % Expand out the local and global indicies for this compressed item
    
    global_row = elementDOF*(elem-1) + [1:elementDOF]';
    
    for d=1:num_dims
        myDeg = pde.deg;
        Index_I{d} = (element_idx1D_D{d}-1)*myDeg + [1:myDeg]';
    end
    
    for j=1:num_connected
              
        if opts.use_connectivity
            connected_col_j = connectivity{elem}(j);  
        else
            connected_col_j = j;
        end   
        
        for d=1:num_dims
            connected_idx1D_D{d} = A_data.element_local_index_D{d}(connected_col_j);
        end
        
        % Expand out the global col indicies for this compressed
        % connected item.
        
        global_col = elementDOF*(connected_col_j-1) + [1:elementDOF]';
        
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
            
            if opts.build_A
                
                %%
                % Apply krond to return A (for hand coded implicit time advance)
                
                view = krond(num_dims,kronMatList);
                A(global_row,global_col) = A(global_row,global_col) + view;

            else
                
                %%
                % Apply kron_mult to return A*Y (explicit time advance)
                X = f(global_col);
                if use_kronmultd
                    Y = kron_multd(num_dims,kronMatList,X);
                else
                    Y = kron_multd_full(num_dims,kronMatList,X);
                end
                
                use_globalRow = 0;
                if (use_globalRow)
                    ftmpA(global_row) = ftmpA(global_row) + Y;
                else
                    % ------------------------------------------------------
                    % globalRow = elementDOF*(workItem-1) + [1:elementDOF]';
                    % ------------------------------------------------------
                    i1 = elementDOF*(elem-1) + 1;
                    i2 = elementDOF*(elem-1) + elementDOF;
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
            
            ALHS(global_row,global_col) = ALHS(global_row,global_col) + krond(num_dims,kronMatListLHS);
            
        end
        
        
        %%
        % Overwrite previous approach with PDE spec approch
        ftmp = ftmpA;
                
    end
    
    assert(elem==elem);
    
end
if opts.build_A && ~opts.quiet; disp('DONE'); end

end