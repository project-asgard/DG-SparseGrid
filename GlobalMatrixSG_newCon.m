function A_encode=GlobalMatrixSG_newCon(A,B,HASH,Lev,Deg)
%Global Matrix construction from Matrices A1,...,ADim
% This code is to generate global matrix
% only works for
% compression = 3
% By Lin, 01/02/2019
%------------------------------------------------------------

global hash_format

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
    % General Dim Method :: use recursive way to generate 
    % needs further work
    Cel1 = repmat(Cell_loc{1}(:),1,numel(Cell_loc{2}));
    Cel2 = repmat(Cell_loc{2}(:),numel(Cell_loc{1}),1);
    Cel1 = Cel1';
    Cel2 = Cel2';
    
    % Compute their index in Hash table
    key = zeros(numel(Cel1),2*Dim);
    key(1:end,1:Dim) = repmat(Lev_loc,numel(Cel1),1);
    key(:,Dim+1:2*Dim) = [Cel1(:),Cel2(:)];
    
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
        end
        
        tmpA=A( Deg*(index_I{1}(1)-1)+1:Deg*(index_I{1}(end)),...
            Deg*(index_J{1}(1)-1)+1:Deg*(index_J{1}(end)) );
        rA = size(tmpA,1);
        cA = size(tmpA,2);
        
        tmpB=B( Deg*(index_I{2}(1)-1)+1:Deg*(index_I{2}(end)),...
            Deg*(index_J{2}(1)-1)+1:Deg*(index_J{2}(end)) );
        rB = size(tmpB,1);
        cB = size(tmpB,2);
        
        % save sub-matries for A and B
        A_encode{count}.A1=tmpA;
        A_encode{count}.A2=tmpB;        
        
        % compute the row index
        tmp = Deg^2*(IndexI(:)-1)+[1:Deg^2];
        tmp = tmp';
        rIndex = kron_split( rA,  Deg*ones(rB/Deg,1) );

        IndexI = tmp(:);
        IndexI(rIndex(:)) = IndexI;
        
        % compute the column index
        tmp = Deg^2*(IndexJ(:)-1)+[1:Deg^2];
        tmp = tmp';
        cIndex = kron_split( cA,  Deg*ones(cB/Deg,1) );

        IndexJ = tmp(:);
        IndexJ(cIndex(:)) = IndexJ;
        
        A_encode{count}.IndexI = IndexI;
        A_encode{count}.IndexJ = IndexJ;
        
        
        count = count+1;
    end
end

