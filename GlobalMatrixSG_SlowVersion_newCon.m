%function [A_encode,A_Data] = GlobalMatrixSG_SlowVersion(coeffMat1,coeffMat2,HASHInv,connectivity,Deg)

function [A_Data] = GlobalMatrixSG_SlowVersion_newCon(HASH,Lev,Deg,compression)

% Global Matrix construction from the coefficient matricies coeffMat1 and
% coeffMat2 by looping over each Lev (i.e., Lev1 + Lev2 < = MaxLev).
% Difference from the Previous Version: Avoid Connectivity
% Only work for compression == 4
%--------------------------------------------------------------------------

global hash_format


if compression ~=4
    error('only can compute compression eq 4. Check value for compression!')
end


Dim = 2;


% All the possible combination of Lev for variables (x1,...xDim)
ComLev = perm_leq( Dim, Lev);
ComLevIndex = [];

for i = 1:size(ComLev,1)
    
    Lev_loc = ComLev(i,:);
    
    for d = 1:Dim
        Cell_loc{d}=[0:2^max(Lev_loc(d)-1,0)-1];
        index = LevCell2index(Lev_loc(d),Cell_loc{d});
        ComLevIndex{i}.lid{d} = sort(index);
    end
    
    % compute all the Cel corresponding to Lev
    % General Dim Method 
    CelMat = CelRepeatGen(Cell_loc);
    
    
    % Compute their index in Hash table
    % for general dim
    key = zeros(size(CelMat,1),2*Dim);
    key(1:end,1:Dim) = repmat(Lev_loc,size(CelMat,1),1);
    key(:,Dim+1:2*Dim) = CelMat;
    
    index = zeros(size(key,1),1);
    
    % find the index from HashTable
    for j =1:size(key,1)
        index(j) = HASH.(sprintf(hash_format,key(j,:)));
    end
    
    ComLevIndex{i}.gid = index;
    clear key Cel1 Cel2 index
end

N = size(ComLevIndex,1);

dofCnt = 1;
% Allocate element arrays

element_global_row_index  = zeros(N,1);
element_local_1_index     = zeros(N,1);
element_local_2_index     = zeros(N,1);

% begin assembling the global matrix
for i = 1:size(ComLev,1)
    for j = 1:size(ComLev,1)
        IndexI = ComLevIndex{i}.gid;
        IndexJ = ComLevIndex{j}.gid;
        
        for d = 1:Dim
            index_I{d} = ComLevIndex{i}.lid{d};
            index_J{d} = ComLevIndex{j}.lid{d};

            rSize{d} = numel(index_I{d});
            cSize{d} = numel(index_J{d});
        end
        
        
        element_idx1D_1 = index_I{1};
        element_idx1D_2 = index_I{2};
        

    
        % !!
        % compute the row index - needs more work for Dim > 2
        tmp = Deg^Dim*(IndexI(:)-1)+[1:Deg^Dim];
        tmp = tmp';
        rIndex = kron_split( Deg*rSize{1},  Deg*ones(rSize{2},1) );
        
        IndexI = tmp(:);
        IndexI(rIndex(:)) = IndexI;
        
        % compute the column index
        tmp = Deg^Dim*(IndexJ(:)-1)+[1:Deg^Dim];
        tmp = tmp';
        cIndex = kron_split( Deg*cSize{1},  Deg*ones(cSize{2},1) );
        
        IndexJ = tmp(:);
        IndexJ(cIndex(:)) = IndexJ;
        
                    % Store the element data in arrays
    element_global_row_index(dofCnt) = IndexI;
    element_local_1_index(dofCnt) = element_idx1D_1;
    element_local_2_index(dofCnt) = element_idx1D_2;
%         
%         A_encode{count}.IndexI = IndexI;
%         A_encode{count}.IndexJ = IndexJ;
        
        dofCnt = dofCnt+1;
    end
end

% Wrap the arrays up into a struct just for ease of passing around.

A_Data.element_global_row_index = element_global_row_index;
A_Data.element_local_1_index = element_local_1_index;
A_Data.element_local_2_index = element_local_2_index;

% Allocate connected element arraysn (these won't be filled, and we remove
% the extra zeros after their construction).

A_Data.connected_global_col_index = connected_global_col_index(1:sum(element_n_connected));
A_Data.connected_local_1_index = connected_local_1_index(1:sum(element_n_connected));
A_Data.connected_local_2_index = connected_local_2_index(1:sum(element_n_connected));



end

function CelMat = CelRepeatGen(CelLoc)
    Dim = numel(CelLoc);
    CelMat = CelLoc{1}(:);
    for d = 1:Dim-1
        tmpCel1 = repmat(CelMat,numel(CelLoc{d+1}),1);
        tmpCel2 = repmat(CelLoc{d+1}(:),1,size(CelMat,1));
        tmpCel2 = tmpCel2';

        CelMat = [tmpCel1,tmpCel2(:)];
    end
end

% count = 1;
% 
% % Use tensor product encoding over Deg.
% 
% dofCnt = 1;
% conCnt = 1;
% 
% 
% % Allocate connected element arraysn (these won't be filled, and we remove
% % the extra zeros after their construction).
% 
% connected_global_col_index  = zeros(N*N,1);
% connected_local_1_index     = zeros(N*N,1);
% connected_local_2_index     = zeros(N*N,1);
% 
% for workItem = 1:N
%     
%     % Get the coordinates in the basis function space for myRow (this
%     % element). (Lev1,Lev2,Cel1,Cel2,idx1D_1,idx1D_2) Lev1,Lev2,Cel1,Cel2
%     % are NOT used here.
%     
%     thisRowBasisCoords = HASHInv{workItem};
%     
%     % Get the 1D indexes into the [lev,pos] space for this element (row)
%     
%     element_idx1D_1 = thisRowBasisCoords(5);
%     element_idx1D_2 = thisRowBasisCoords(6);
%     
%     % Get the global index of non-zero (connected) columns for this row
%     
%     connectedCols = connectivity{workItem};
%     
%     % Get the local (basis) coords the connected elements
%     %
%     % Hash :    local coords  -> global index HashInv:  global index  ->
%     % local coords
%     
%     connectedColsBasisCoords = [HASHInv{connectedCols}];
%     
%     % Get the 1D indices into the [lev,pos] space for the connected
%     % elements (cols)
%     
%     connected_idx1D_1 = connectedColsBasisCoords(5:6:end);
%     connected_idx1D_2 = connectedColsBasisCoords(6:6:end);
%     
%     % Store the element data in arrays
%     element_global_row_index(dofCnt) = workItem;
%     element_local_1_index(dofCnt) = element_idx1D_1;
%     element_local_2_index(dofCnt) = element_idx1D_2;
%     
%     nConnections = 0;
%     % Loop over connected elements
%     for jjj = 1:size(connected_idx1D_1,2)
%         
%         % Store the connected data in arrays
%         connected_global_col_index(conCnt) = connectedCols(jjj);
%         connected_local_1_index(conCnt) = connected_idx1D_1(jjj);
%         connected_local_2_index(conCnt) = connected_idx1D_2(jjj);
%         
%         conCnt = conCnt+1;
%         nConnections = nConnections + 1;
%         
%     end
%     
%     element_n_connected(dofCnt) = nConnections;
%     
%     dofCnt = dofCnt + 1;
%     
% end
% 
% % Wrap the arrays up into a struct just for ease of passing around.
% 
% A_Data.element_global_row_index = element_global_row_index;
% A_Data.element_local_1_index = element_local_1_index;
% A_Data.element_local_2_index = element_local_2_index;
% A_Data.element_n_connected = element_n_connected;
% 
% % Allocate connected element arraysn (these won't be filled, and we remove
% % the extra zeros after their construction).
% 
% A_Data.connected_global_col_index = connected_global_col_index(1:sum(element_n_connected));
% A_Data.connected_local_1_index = connected_local_1_index(1:sum(element_n_connected));
% A_Data.connected_local_2_index = connected_local_2_index(1:sum(element_n_connected));
% 
% 
% 
% end
