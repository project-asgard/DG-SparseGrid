function [M] = moment_reduced_matrix(opts,pde,A_data,hash_table,moment_idx)
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CURRENTLY ONLY WORKS FOR 1x-1v,
%         1x-2v, or 1x-3v system
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% This function produces a matrix M that is the transformation matrix that
% maps sparse-grid functions in (x,v) to DG functions only in x.  This
% transformation is given by 
%          f_x(x) = \int_v f(x,v)*g(v) dv
% where g(v) is a moment function with index moment_idx
% 
% In 1x-1v systems, the output function is full grid, but this will break
% down in multiple dimensions.   

num_dims = numel(pde.dimensions);
deg = opts.deg;
num_ele = numel(A_data.element_local_index_D{1});

if num_dims == 2 % 1x-1v

    x_dim = 1; %Hardcoded for now.  Needs to change in multi-dim
    v_dim = 2; %Hardcoded for now.  Needs to change in multi-dim
    g_vec = pde.moments{moment_idx}.fList{v_dim};

    %Get row-dim of M -- Assuming a full grid output here
    rows = deg*2^ceil(log2(max(A_data.element_local_index_D{x_dim})));

    %Containers for sparse matrix
    I = zeros(deg^2*num_ele,1); %Square is hardcoded.  Needs to change 
    J = zeros(deg^2*num_ele,1);
    S = zeros(deg^2*num_ele,1);

    count = 1;
    k = 1:deg; %Not going to loop over degree
    for i=1:num_ele
        for j=1:deg
            temp = ((A_data.element_local_index_D{x_dim}(i)-1)*deg + j) + 0*k;
            I(count:count+deg-1) = temp';

            temp = (i-1)*deg^2 + (j-1)*deg + k;
            J(count:count+deg-1) = temp';

            temp = g_vec((A_data.element_local_index_D{v_dim}(i)-1)*deg + k);
            S(count:count+deg-1) = temp';

            count = count + deg;
        end
    end

    M = sparse(I,J,S,rows,deg^2*num_ele);

elseif num_dims == 3 % 1x-2v
    
    x_dim   = 1; %Hardcoded for now.  Needs to change in multi-dim
    v_dim_1 = 2; %Hardcoded for now.  Needs to change in multi-dim
    v_dim_2 = 3; %Hardcoded for now.  Needs to change in multi-dim
    
    %Get row-dim of M -- Assuming a full grid output here
    rows = deg*2^ceil(log2(max(A_data.element_local_index_D{x_dim})));

    %Containers for sparse matrix
    I = zeros(deg^3*num_ele,1); %Square is hardcoded.  Needs to change 
    J = zeros(deg^3*num_ele,1);
    S = zeros(deg^3*num_ele,1);
    
    g_vec_1 = pde.moments{moment_idx}.fList{v_dim_1};
    g_vec_2 = pde.moments{moment_idx}.fList{v_dim_2};
    
    count = 0;
    for i=1:num_ele %loop through all elements lit up
        
        %Get local dof for element
        l_dof_x = A_data.element_local_index_D{x_dim  }(i)-1;
        l_dof_v = [A_data.element_local_index_D{v_dim_1}(i)-1,...
                   A_data.element_local_index_D{v_dim_2}(i)-1];
                 
        for xdeg=1:deg
            %The row of the matrix is given by the x-dof
            row_idx = l_dof_x*deg+xdeg;
            for vdeg1=1:deg
                for vdeg2=1:deg
                    %Get column entry based on degree index
                    col_idx = (i-1)*deg^3 + (xdeg-1)*deg^2 + (vdeg1-1)*deg + vdeg2;
                    %Get entry (which is the product of the 1d elements of
                    %g_vec)
                    val = g_vec_1(l_dof_v(1)*deg+vdeg1) * ...
                          g_vec_2(l_dof_v(2)*deg+vdeg2);
                    I(count+1) = row_idx;
                    J(count+1) = col_idx;
                    S(count+1) = val;
                    count = count + 1;
                end
            end
        end
    end
    
    %Construct sparse matrix
    M = sparse(I,J,S,rows,deg^3*num_ele);
   

elseif num_dims == 4 % 1x-3v

    x_dim   = 1; %Hardcoded for now.  Needs to change in multi-dim
    v_dim_1 = 2; %Hardcoded for now.  Needs to change in multi-dim
    v_dim_2 = 3; %Hardcoded for now.  Needs to change in multi-dim
    v_dim_3 = 4; %Hardcoded for now.  Needs to change in multi-dim
    
    %Get row-dim of M -- Assuming a full grid output here
    rows = deg*2^ceil(log2(max(A_data.element_local_index_D{x_dim})));

    %Containers for sparse matrix
    I = zeros(deg^4*num_ele,1);
    J = zeros(deg^4*num_ele,1);
    S = zeros(deg^4*num_ele,1);
    
    g_vec_1 = pde.moments{moment_idx}.fList{v_dim_1};
    g_vec_2 = pde.moments{moment_idx}.fList{v_dim_2};
    g_vec_3 = pde.moments{moment_idx}.fList{v_dim_3};
    
    count = 0;
    for i=1:num_ele %loop through all elements lit up
        
        %Get local dof for element
        l_dof_x = A_data.element_local_index_D{x_dim  }(i)-1;
        l_dof_v = [A_data.element_local_index_D{v_dim_1}(i)-1,...
                   A_data.element_local_index_D{v_dim_2}(i)-1,...
                   A_data.element_local_index_D{v_dim_3}(i)-1];
                 
        for xdeg=1:deg
            %The row of the matrix is given by the x-dof
            row_idx = l_dof_x*deg+xdeg;
            for vdeg1=1:deg
                for vdeg2=1:deg
                    for vdeg3=1:deg
                        %Get column entry based on degree index
                        col_idx = (i-1)*deg^4 + (xdeg-1)*deg^3 + (vdeg1-1)*deg^2 + (vdeg2-1)*deg + vdeg3;
                        %Get entry (which is the product of the 1d elements of
                        %g_vec)
                        val = g_vec_1(l_dof_v(1)*deg+vdeg1) * ...
                              g_vec_2(l_dof_v(2)*deg+vdeg2) * ...
                              g_vec_3(l_dof_v(3)*deg+vdeg3);
                        I(count+1) = row_idx;
                        J(count+1) = col_idx;
                        S(count+1) = val;
                        count = count + 1;
                    end
                end
            end
        end
    end
    
    %Construct sparse matrix
    M = sparse(I,J,S,rows,deg^4*num_ele);

end


end

