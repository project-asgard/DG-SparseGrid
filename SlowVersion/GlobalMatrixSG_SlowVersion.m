%function [A_encode,A_Data] = GlobalMatrixSG_SlowVersion(coeffMat1,coeffMat2,HASHInv,connectivity,Deg)

function [A_Data] = GlobalMatrixSG_SlowVersion(pde,HASHInv,connectivity,Deg,compression)

% Global Matrix construction from the coefficient matricies coeffMat1 and
% coeffMat2 by looping over each grid point (i.e., elements of the hash).
% Each grid point is represented by a row in the global system matrix A,
% with the non-zero elements of that row being the other grid points to
% which the point is connected.

N = size(HASHInv,2);

nDims = numel(pde.dimensions);

if compression < 3
    
    % No tensor product encoding of A.
    
    nDOF = N * Deg^2; % For a 2D problem. Deg^3 for 3D etc.
    
    dofCnt = 1;
    conCnt = 1;
    
    % Allocate element arrays
    
    element_global_row_index  = zeros(nDOF,1);
    element_local_1_index     = zeros(nDOF,1);
    element_local_2_index     = zeros(nDOF,1);
    element_n_connected       = zeros(nDOF,1);
    
    % Allocate connected element arraysn (these won't be filled, and we remove
    % the extra zeros after their construction).
    
    connected_global_col_index  = zeros(nDOF*nDOF,1);
    connected_local_1_index     = zeros(nDOF*nDOF,1);
    connected_local_2_index     = zeros(nDOF*nDOF,1);
    
    for workItem = 1:N
        
        % Get the coordinates in the basis function space for myRow (this
        % element). (Lev1,Lev2,Cel1,Cel2,idx1D_1,idx1D_2) Lev1,Lev2,Cel1,Cel2
        % are NOT used here.
        
        thisRowBasisCoords = HASHInv{workItem};
        
        % Get the 1D indexes into the [lev,pos] space for this element (row)
        
        element_idx1D_1 = thisRowBasisCoords(5);
        element_idx1D_2 = thisRowBasisCoords(6);
        
        % Get the global index of non-zero (connected) columns for this row
        
        connectedCols = connectivity{workItem};
        
        % Get the local (basis) coords the connected elements
        %
        % Hash :    local coords  -> global index HashInv:  global index  ->
        % local coords
        
        connectedColsBasisCoords = [HASHInv{connectedCols}];
        
        % Get the 1D indices into the [lev,pos] space for the connected
        % elements (cols)
        
        connected_idx1D_1 = connectedColsBasisCoords(5:6:end);
        connected_idx1D_2 = connectedColsBasisCoords(6:6:end);
        
        % Loop over dim1 Deg for this element
        for k1 = 1:Deg
            
            
            % Get dim1 1D index into [lev,pos,deg] space (for this element)
            index_I1 = (element_idx1D_1-1)*Deg+k1;
            
            
            % Loop over dim2 Deg for this element
            for k2 = 1:Deg
                
                
                % Get dim2 1D index into [lev,pos,deg] space (for this element)
                index_I2 = (element_idx1D_2-1)*Deg+k2;
                
                % Get the global row index for this element
                globalRow = Deg^2*(workItem-1)+Deg*(k1-1)+k2;
                
                % Store the element data in arrays
                element_global_row_index(dofCnt) = globalRow;
                element_local_1_index(dofCnt) = index_I1;
                element_local_2_index(dofCnt) = index_I2;
                
                nConnections = 0;
                % Loop over connected elements
                for jjj = 1:size(connected_idx1D_1,2)
                    
                    % Loop over dim1 Deg for the connected elements
                    for kk1 = 1:Deg
                        
                        
                        % Get dim1 1D index into [lev,pos,deg] space (for
                        % connected element)
                        index_J1 = (connected_idx1D_1(jjj)-1)*Deg+kk1;
                        
                        
                        % Loop over dim2 Deg for connected elements
                        for kk2 = 1:Deg
                            
                            
                            % Get dim2 1D index into [lev,pos,deg] space (for
                            % connected element)
                            index_J2 = (connected_idx1D_2(jjj)-1)*Deg+kk2;
                            
                            % Get the global column index for this connected
                            % element.
                            globalCol = Deg^2*(connectedCols(jjj)-1)+Deg*(kk1-1)+kk2;
                            
                            % Store the connected data in arrays
                            connected_global_col_index(conCnt) = globalCol;
                            connected_local_1_index(conCnt) = index_J1;
                            connected_local_2_index(conCnt) = index_J2;
                            
                            %                         tmpA = coeffMat1(index_I1,index_J1);
                            %                         tmpB = coeffMat2(index_I2,index_J2);
                            %
                            %                         A_encode{conCnt}.IndexI = globalRow;
                            %
                            %                         A_encode{conCnt}.A1=tmpA;
                            %                         A_encode{conCnt}.A2=tmpB;
                            %
                            %
                            %                         A_encode{conCnt}.IndexJ = globalCol;
                            
                            conCnt = conCnt+1;
                            nConnections = nConnections + 1;
                        end
                    end
                end
                
                element_n_connected(dofCnt) = nConnections;
                
                dofCnt = dofCnt + 1;
                
            end
        end
    end
    
    
    % Wrap the arrays up into a struct just for ease of passing around.
    
    A_Data.element_global_row_index = element_global_row_index;
    A_Data.element_local_1_index = element_local_1_index;
    A_Data.element_local_2_index = element_local_2_index;
    A_Data.element_n_connected = element_n_connected;
    
    % Allocate connected element arraysn (these won't be filled, and we remove
    % the extra zeros after their construction).
    
    A_Data.connected_global_col_index = connected_global_col_index(1:sum(element_n_connected));
    A_Data.connected_local_1_index = connected_local_1_index(1:sum(element_n_connected));
    A_Data.connected_local_2_index = connected_local_2_index(1:sum(element_n_connected));
    
    
elseif compression == 4
    
    % Use tensor product encoding over Deg.
    
    dofCnt = 1;
    conCnt = 1;
    
    % Allocate element arrays
    
    element_global_row_index  = zeros(N,1);
    element_local_1_index     = zeros(N,1);
    element_local_2_index     = zeros(N,1);
    for d=1:nDims
        element_local_index_D{d} = zeros(N,1);
    end
    element_n_connected       = zeros(N,1);
    
    % Allocate connected element arraysn (these won't be filled, and we remove
    % the extra zeros after their construction).
    
    connected_global_col_index  = zeros(N*N,1);
    connected_local_1_index     = zeros(N*N,1);
    connected_local_2_index     = zeros(N*N,1);
    for d=1:nDims
        connected_local_index_D{d} = zeros(N*N,1);
    end
    
    for workItem = 1:N
        
        % Get the coordinates in the basis function space for myRow (this
        % element). (Lev1,Lev2,Cel1,Cel2,idx1D_1,idx1D_2) Lev1,Lev2,Cel1,Cel2
        % are NOT used here.
        
        thisRowBasisCoords = HASHInv{workItem};
        
        % Get the 1D indexes into the [lev,pos] space for this element (row)
        
        element_idx1D_1 = thisRowBasisCoords(5);
        element_idx1D_2 = thisRowBasisCoords(6);
        for d=1:nDims
            element_idx1D_D{d} = thisRowBasisCoords(nDims*2+d);
        end
        
        % Get the global index of non-zero (connected) columns for this row
        
        connectedCols = connectivity{workItem};
        
        % Get the local (basis) coords the connected elements
        %
        % Hash : local coords  -> global index 
        % HashInv:  global index -> local coords
        
        connectedColsBasisCoords = [HASHInv{connectedCols}];
        
        % Get the 1D indices into the [lev,pos] space for the connected
        % elements (cols)
        
        connected_idx1D_1 = connectedColsBasisCoords(5:6:end);
        connected_idx1D_2 = connectedColsBasisCoords(6:6:end);
        for d=1:nDims    
            % Recall ...
            % 1D : (lev1,cell1,idx1D1)
            % 2D : (lev1,lev2,cell1,cell2,idx1D1,idx1D2)
            % 3D : (lev1,lev2,lev3,cell1,cell2,cell3,idx1D1,idx1D2,idx1D3)
            % such that we have the following indexing generalized to
            % dimension ...
            connected_idx1D_D{d} = connectedColsBasisCoords(2*nDims+d:3*nDims:end);
        end
        
        % Store the element data in arrays
        element_global_row_index(dofCnt) = workItem;
        
        element_local_1_index(dofCnt) = element_idx1D_1;
        element_local_2_index(dofCnt) = element_idx1D_2;
        for d=1:nDims       
            element_local_index_D{d}(dofCnt) = element_idx1D_D{d};
        end
        
        nConnections = 0;
        % Loop over connected elements
        assert(size(connected_idx1D_D{1},2) == numel(connectedCols));
%         for jjj = 1:size(connected_idx1D_1,2)
        for jjj = 1:numel(connectedCols)
            
            % Store the connected data in arrays
            connected_global_col_index(conCnt) = connectedCols(jjj);
            
            connected_local_1_index(conCnt) = connected_idx1D_1(jjj);
            connected_local_2_index(conCnt) = connected_idx1D_2(jjj);
            for d=1:nDims
                connected_local_index_D{d}(conCnt) = connected_idx1D_D{d}(jjj);
            end
            
            conCnt = conCnt+1;
            nConnections = nConnections + 1;

        end
        
        element_n_connected(dofCnt) = nConnections;
        
        dofCnt = dofCnt + 1;
   
    end
    
    % Wrap the arrays up into a struct just for ease of passing around.
    
    A_Data.element_global_row_index = element_global_row_index;
    A_Data.element_local_1_index = element_local_1_index;
    A_Data.element_local_2_index = element_local_2_index;
    A_Data.element_local_index_D = element_local_index_D;
    A_Data.element_n_connected = element_n_connected;
    
    % Allocate connected element arraysn (these won't be filled, and we remove
    % the extra zeros after their construction).
    
    A_Data.connected_global_col_index = connected_global_col_index(1:sum(element_n_connected));
    A_Data.connected_local_1_index = connected_local_1_index(1:sum(element_n_connected));
    A_Data.connected_local_2_index = connected_local_2_index(1:sum(element_n_connected));
    for d=1:nDims
        A_Data.connected_local_index_D{d} = connected_local_index_D{d}(1:sum(element_n_connected));
    end
    
end

end
