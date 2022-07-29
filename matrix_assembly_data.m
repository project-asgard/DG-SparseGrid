function [ A_Data ] = matrix_assembly_data( unknown )

    % Global Matrix construction from the coefficient matricies coeffMat1 and
    % coeffMat2 by looping over each grid point (i.e., elements of the hash).
    % Each grid point is represented by a row in the global system matrix A,
    % with the non-zero elements of that row being the other grid points to
    % which the point is connected.

    Ne = numel(unknown.hash_table.elements_idx);
    N = Ne;

    nDims = numel(unknown.dimensions);

    % Use tensor product encoding over deg.

    dofCnt = 1;
    conCnt = 1;

    % Allocate element arrays

    element_global_row_index = zeros(N,1);

    for d=1:nDims
        element_local_index_D{d} = zeros(N,1);
    end

    for workItem = 1 : N
    
        % Get the coordinates in the basis function space for myRow (this
        % element). (Lev1,Lev2,Cel1,Cel2,idx1D_1,idx1D_2) Lev1,Lev2,Cel1,Cel2
        % are NOT used here.
    
        % Get the 1D indexes into the [lev,pos] space for this element (row)
    
        for d=1 : nDims
            
            IDlev = unknown.hash_table.elements.lev_p1(unknown.hash_table.elements_idx(workItem),d)-1;
            IDpos = unknown.hash_table.elements.pos_p1(unknown.hash_table.elements_idx(workItem),d)-1;
            IDe = lev_cell_to_1D_index(IDlev,IDpos);
            element_idx1D_D{d} = IDe;
            
        end
    
        % Store the element data in arrays
        element_global_row_index(dofCnt) = workItem;
    
        for d = 1 : nDims
            
            element_local_index_D{d}(dofCnt) = element_idx1D_D{d};
            
        end
        
        dofCnt = dofCnt + 1;
        
    end

    % Wrap the arrays up into a struct just for ease of passing around.

    A_Data.element_global_row_index = element_global_row_index;
    A_Data.element_local_index_D    = element_local_index_D;

end
