function A_encode=GlobalMatrixSG_SlowVersion(coeffMat1,coeffMat2,HASHInv,connectivity,Deg)
% Global Matrix construction from the coefficient matricies coeffMat1 and
% coeffMat2 by looping over each grid point (i.e., elements of the hash).
% Each grid point is represented by a row in the global system matrix A,
% with the non-zero elements of that row being the other grid points to
% which the point is connected.

nRows = size(HASHInv,2);

count = 1;

for thisRow = 1:nRows
    
    % Get the coordinates in the basis function space for myRow (this
    % element). (Lev1,Lev2,Cel1,Cel2,idx1D_1,idx1D_2)
    % Lev1,Lev2,Cel1,Cel2 are NOT used here. 
    
    thisRowBasisCoords = HASHInv{thisRow};
    
    % Get the 1D indexes into the [lev,pos] space for this element (row)
    
    idx1D_1 = thisRowBasisCoords(5);
    idx1D_2 = thisRowBasisCoords(6);
    
    % Get the global index of non-zero (connected) columns for this row
    
    connectedCols = connectivity{thisRow};
    
    % Get the local (basis) coords the connected elements
    %
    % Hash :    local coords  -> global index
    % HashInv:  global index  -> local coords
    
    connectedColsBasisCoords = [HASHInv{connectedCols}];
    
    % Get the 1D indices into the [lev,pos] space for the connected elements (cols)
    
    connected_idx1D_1 = connectedColsBasisCoords(5:6:end);
    connected_idx1D_2 = connectedColsBasisCoords(6:6:end);
    
    % Loop over dim1 Deg for this element
    for k1 = 1:Deg
        
        
        % Get dim1 1D index into [lev,pos,deg] space (for this element)
        index_I1 = (idx1D_1-1)*Deg+k1;
        
        
        % Loop over dim2 Deg for this element
        for k2 = 1:Deg
            
            
            % Get dim2 1D index into [lev,pos,deg] space (for this element)
            index_I2 = (idx1D_2-1)*Deg+k2;
            
            % Get the global row index for this element
            globalRow = Deg^2*(thisRow-1)+Deg*(k1-1)+k2;
            
            
            % Loop over connected elements
            for jjj = 1:size(connected_idx1D_1,2)
                
                % Loop over dim1 Deg for the connected elements
                for kk1 = 1:Deg
                    
                    
                    % Get dim1 1D index into [lev,pos,deg] space 
                    % (for connected element)
                    index_J1 = (connected_idx1D_1(jjj)-1)*Deg+kk1;
                    
                    
                    % Loop over dim2 Deg for connected elements
                    for kk2 = 1:Deg
                       
                        
                        % Get dim2 1D index into [lev,pos,deg] space
                        % (for connected element)
                        index_J2 = (connected_idx1D_2(jjj)-1)*Deg+kk2;
                        
                        
                        % Get the global column index for this connected
                        % element.
                        globalCol = Deg^2*(connectedCols(jjj)-1)+Deg*(kk1-1)+kk2;
                        
                        tmpA = coeffMat1(index_I1,index_J1);
                        tmpB = coeffMat2(index_I2,index_J2);
                        
                        A_encode{count}.IndexI = globalRow;
                        
                        A_encode{count}.A1=tmpA;
                        A_encode{count}.A2=tmpB;
                        
                        
                        A_encode{count}.IndexJ = globalCol;
                        
                        count = count+1;
                    end
                end
                
            end
            
        end
    end
    
     
end
