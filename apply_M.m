function [ftmp] = apply_M (pde,opts,A_data,f)

%-----------------------------------
% Multiply Matrix M by a vector f
%-----------------------------------
dof = size(f,1);
use_sparse_ftmp = 0;
if (use_sparse_ftmp)
    ftmp=sparse(dof,1);
else
    ftmp = zeros(dof,1);
end
use_kronmultd = 1;

deg = opts.deg;

num_dims      = numel(pde.dimensions);

%%
% Tensor product encoding over DOF within an element, i.e., over "deg" (A_Data),
% i.e., tmpA and tmpB are deg_1 x deg_2 x deg_D matricies

num_elem = numel(A_data.element_global_row_index);

ftmpA = ftmp;

element_DOF = deg^num_dims;

total_DOF = num_elem * element_DOF;

if opts.use_connectivity
    connectivity = pde.connectivity;
    num_A = 0;
    for i=1:num_elem
        num_A = num_A + numel(pde.connectivity{i});
    end
    num_A = num_A * element_DOF^2;
else
    num_A = total_DOF * total_DOF; 
end

cnt = 1;

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
    
    global_row = element_DOF*(elem-1) + [1:element_DOF]';
%     global_1D_row = deg*(elem-1) + [1:deg]';
    
    for d=1:num_dims
        myDeg = opts.deg;
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
        
        global_col = element_DOF*(connected_col_j-1) + [1:element_DOF]';
        
        for d=1:num_dims
            myDeg = opts.deg;
            Index_J{d} = (connected_idx1D_D{d}-1)*myDeg + [1:myDeg]';
        end
        
        %%
        % Apply operator matrices to present state using the pde spec
        % Y = A * X
        % where A is tensor product encoded.
 
        if opts.build_A && opts.use_sparse_A
            num_view = element_DOF * element_DOF;
            [gr,gc] = meshgrid(global_col,global_row);
            A_s1(cnt:cnt+num_view-1) = gr(:);
            A_s2(cnt:cnt+num_view-1) = gc(:);
        end
        
        %%
        % Construct the list of matrices for the kron_mult for this
        % operator (which has dimension many entries).
        for d=1:num_dims
            idx_i = Index_I{d};
            idx_j = Index_J{d};
            tmp = pde.dimensions{d}.mass_mat;
            kron_mat_list{d} = tmp(idx_i,idx_j); % List of tmpA, tmpB, ... tmpD used in kron_mult
        end


        %%
        % Apply kron_mult to return A*Y (explicit time advance)
        X = f(global_col);
        if use_kronmultd
            Y = kron_multd(num_dims,kron_mat_list,X);
        else
            Y = kron_multd_full(num_dims,kron_mat_list,X);
        end

        use_globalRow = 0;
        if (use_globalRow)
            ftmpA(global_row) = ftmpA(global_row) + Y;
        else
            % ------------------------------------------------------
            % globalRow = elementDOF*(workItem-1) + [1:elementDOF]';
            % ------------------------------------------------------
            i1 = element_DOF*(elem-1) + 1;
            i2 = element_DOF*(elem-1) + element_DOF;
            ftmpA(i1:i2) = ftmpA(i1:i2) + Y;
        end
        
        %%
        % Overwrite previous approach with PDE spec approch
        ftmp = ftmpA;
        
        if opts.use_sparse_A; cnt = cnt + num_view; end
                
    end
    
    assert(elem==elem);
    
end

end
