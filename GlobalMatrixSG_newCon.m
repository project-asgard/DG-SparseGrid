function A_encode=GlobalMatrixSG_newCon(A,B,HASH,Lev,Deg,gridType)
%Global Matrix construction from Matrices A1,...,ADim
% This code is to generate global matrix
% only works for
% compression = 3
% Difference with Previous Version: Avoid the connectivity
% Needs more work for Line 76 and 84
%--------------------------------------------------------------------------

global hash_format

Dim = 2;
Matrix{1} = A;
Matrix{2} = B;

is_sparse_grid = strcmp( gridType, 'SG');

% All the possible combination of Lev for variables (x1,...xDim)
if (is_sparse_grid),
    ComLev = perm_leq( Dim, Lev);
else
    ComLev = perm_max( Dim, Lev);
end
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
count = 1;

% begin assembling the global matrix
for i = 1:size(ComLev,1)
    for j = 1:size(ComLev,1)
        IndexI = ComLevIndex{i}.gid;
        IndexJ = ComLevIndex{j}.gid;
        
        for d = 1:Dim
            index_I{d} = ComLevIndex{i}.lid{d};
            index_J{d} = ComLevIndex{j}.lid{d};
            
            tmpA{d}=Matrix{d}( Deg*(index_I{d}(1)-1)+1:Deg*(index_I{d}(end)),...
                Deg*(index_J{d}(1)-1)+1:Deg*(index_J{d}(end)) );
            rSize{d} = size(tmpA{d},1)/Deg;
            cSize{d} = size(tmpA{d},2)/Deg;
            
            % save sub-matries
            A_encode{count}.A{d} = tmpA{d};
        end
        
        if Dim == 2
            A_encode{count}.A1=A_encode{count}.A{1};
            A_encode{count}.A2=A_encode{count}.A{2};
        end
        
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

        
        A_encode{count}.IndexI = IndexI;
        A_encode{count}.IndexJ = IndexJ;
        
        count = count+1;
    end
end
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

