function [ Af ] = apply_A_term( opts, md_term, input_unknown, output_unknown, input_coeffs )

deg = opts.deg;

%----------------------------------------------
% Multiply Matrix A by Vector f (input_coeffs)
%----------------------------------------------
dof_in =  input_unknown.size();
dof_ot = output_unknown.size();

use_kronmultd = 1;

num_dims_in = numel( input_unknown.dimensions);
num_dims_ot = numel(output_unknown.dimensions);

use_connectivity = 0; % TODO: Connectivity to be added later

%%
% Tensor product encoding over DOF within an element, i.e., over "deg" (A_Data),
% i.e., tmpA and tmpB are deg_1 x deg_2 x deg_D matricies

num_elem_in = numel( input_unknown.A_data.element_global_row_index);
num_elem_ot = numel(output_unknown.A_data.element_global_row_index);

element_DOF_in = deg^num_dims_in;
element_DOF_ot = deg^num_dims_ot;

num_A = dof_in * dof_ot;

Af = zeros(dof_ot,1);

if opts.fast_FG_matrix_assembly
    
    assert( false, 'apply_A_term not implemented for fast_FG_matrix_assembly' )

    opts.use_sparse_A = false;
    
    if num_dims_in == 1
        
        assert(opts.adapt == false,"'fast_FG_matrix_assembly' is not supported for adaptivity and dim = 1");
        %In this case FG == SG
        
        iperm = zeros(dof_in,1);
        perm = zeros(dof_in,1);
        for i=1:deg
            iperm(i:deg:end) = deg*(A_data.element_local_index_D{1}-1)+i;
        end
        perm(iperm) = 1:dof_in;
        
        f_F = input_coeffs(perm);
        Af_F = md_term.terms_1D{1}.mat*f_F;
        Af = Af_F(iperm);
        
        
    elseif num_dims_in == 2
    
        % Get SG <-> FG conversion
        [perm,iperm,pvec] = sg_to_fg_mapping_2d(md_term,opts,A_data);

        f_F = zeros(size(perm,1),1);
        f_F(pvec) = input_coeffs(perm(pvec));
        Af_F = md_term.terms_1D{2}.mat * reshape(f_F,size(md_term.terms_1D{2}.mat,1),[]) * (md_term.terms_1D{1}.mat');
        Af_F = reshape(Af_F,[],1);

        Af = Af_F(iperm);
    
    else
       
       error("'fast_FG_matrix_assembly' only applies to 1 and 2 dimension problems."); 
        
    end
    
else % do not use fast_FG_matrix_assembly
    
    in_element_idx1D_D = cell(num_dims_in,1);
    ot_element_idx1D_D = cell(num_dims_ot,1);

    Index_I = cell(num_dims_in,1);
    Index_J = cell(num_dims_ot,1);

    for elem = 1 : num_elem_ot
        
        for d = 1 : num_dims_ot
            ot_element_idx1D_D{d} = output_unknown.A_data.element_local_index_D{d}(elem);
        end
        
        % Expand out the local and global indicies for this compressed item
        
        global_row = element_DOF_ot*(elem-1) + (1:element_DOF_ot)';
        
        for d = 1 : num_dims_in
            if( d<=num_dims_ot )
                Index_I{d} = (ot_element_idx1D_D{d}-1)*deg + (1:deg)';
            else % TODO: Specialized for (x,v) -> (x).  Needs to be generalized!
                Index_I{d} = 1;
            end
        end
        
        for j = 1 : num_elem_in
            
            for d = 1 : num_dims_in
                in_element_idx1D_D{d} = input_unknown.A_data.element_local_index_D{d}(j);
            end
            
            % Expand out the global col indicies for this compressed
            % connected item.
            
            global_col = element_DOF_in*(j-1) + (1:element_DOF_in)';
            
            for d = 1 : num_dims_in
                Index_J{d} = (in_element_idx1D_D{d}-1)*deg + (1:deg)';
            end
            
            %%
            % Apply operator matrices to present state using the pde spec
            % Y = A * X
            % where A is tensor product encoded.
            
            %%
            % Construct the list of matrices for the kron_mult for this
            % operator (which has dimension many entries).
            for d = 1 : num_dims_in
                idx_i = Index_I{d};
                idx_j = Index_J{d};
                tmp = md_term.terms_1D{d}.mat;
                kron_mat_list{d} = tmp(idx_i,idx_j); % List of tmpA, tmpB, ... tmpD used in kron_mult
            end

            %%
            % Apply kron_mult to return A*Y (explicit time advance)
            X = input_coeffs(global_col);
            if use_kronmultd
                Y = kron_multd(num_dims_in,kron_mat_list,X);
            else
                Y = kron_multd_full(num_dims_in,kron_mat_list,X);
            end

            use_globalRow = 0;
            if (use_globalRow)
                Af(global_row) = Af(global_row) + Y;
            else
                % ------------------------------------------------------
                % globalRow = elementDOF*(workItem-1) + [1:elementDOF]';
                % ------------------------------------------------------
                i1 = element_DOF_ot*(elem-1) + 1;
                i2 = element_DOF_ot*(elem-1) + element_DOF_ot;
                Af(i1:i2) = Af(i1:i2) + Y;
            end
            
        end
        
    end
end

end
