function [ Af ] = apply_A_term( opts, md_term, A_data, f, deg )

%-----------------------------------
% Multiply Matrix A by Vector f
%-----------------------------------
dof = size(f,1);
use_sparse_Af = 0;
if (use_sparse_Af)
    Af=sparse(dof,1);
else
    Af = zeros(dof,1);
end
use_kronmultd = 1;

num_dims = numel(md_term.terms_1D);

use_connectivity = 0; % TODO: Connectivity to be added later

%%
% Tensor product encoding over DOF within an element, i.e., over "deg" (A_Data),
% i.e., tmpA and tmpB are deg_1 x deg_2 x deg_D matricies

num_elem = numel(A_data.element_global_row_index);

element_DOF = deg^num_dims;

total_DOF = num_elem * element_DOF;

if use_connectivity
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

if opts.fast_FG_matrix_assembly
    
    opts.use_sparse_A = false;
    
    if num_dims == 1
        
        assert(opts.adapt == false,"'fast_FG_matrix_assembly' is not supported for adaptivity and dim = 1");
        %In this case FG == SG
        
        iperm = zeros(dof,1);
        perm = zeros(dof,1);
        for i=1:deg
            iperm(i:deg:end) = deg*(A_data.element_local_index_D{1}-1)+i;
        end
        perm(iperm) = 1:dof;
        
        f_F = f(perm);
        Af_F = md_term.terms_1D{1}.mat*f_F;
        Af = Af_F(iperm);
        
        
    elseif num_dims == 2
    
        % Get SG <-> FG conversion
        [perm,iperm,pvec] = sg_to_fg_mapping_2d(md_term,opts,A_data);

        f_F = zeros(size(perm,1),1);
        f_F(pvec) = f(perm(pvec));
        Af_F = md_term.terms_1D{2}.mat * reshape(f_F,size(md_term.terms_1D{1}.mat,1),[]) * (md_term.terms_1D{1}.mat');
        Af_F = reshape(Af_F,[],1);

        Af = Af_F(iperm);
    
    else
       
       error("'fast_FG_matrix_assembly' only applies to 1 and 2 dimension problems."); 
        
    end
    
else % do not use fast_FG_matrix_assembly
    
    for elem=1:num_elem
        
        if use_connectivity
            num_connected = numel(connectivity{elem});
        else
            num_connected = num_elem; % Simply assume all are connected.
        end
        
        for d=1:num_dims
            element_idx1D_D{d} = A_data.element_local_index_D{d}(elem);
        end
        
        % Expand out the local and global indicies for this compressed item
        
        global_row = element_DOF*(elem-1) + [1:element_DOF]';
        
        for d=1:num_dims
            myDeg = deg;
            Index_I{d} = (element_idx1D_D{d}-1)*myDeg + [1:myDeg]';
        end
        
        for j=1:num_connected
            
            if use_connectivity
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
                myDeg = deg;
                Index_J{d} = (connected_idx1D_D{d}-1)*myDeg + [1:myDeg]';
            end
            
            %%
            % Apply operator matrices to present state using the pde spec
            % Y = A * X
            % where A is tensor product encoded.
            
            %%
            % Construct the list of matrices for the kron_mult for this
            % operator (which has dimension many entries).
            for d=1:num_dims
                idx_i = Index_I{d};
                idx_j = Index_J{d};
                tmp = md_term.terms_1D{d}.mat;
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
                Af(global_row) = Af(global_row) + Y;
            else
                % ------------------------------------------------------
                % globalRow = elementDOF*(workItem-1) + [1:elementDOF]';
                % ------------------------------------------------------
                i1 = element_DOF*(elem-1) + 1;
                i2 = element_DOF*(elem-1) + element_DOF;
                Af(i1:i2) = Af(i1:i2) + Y;
            end
            
        end
        
    end
end

end
