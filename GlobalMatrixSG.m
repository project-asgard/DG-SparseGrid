function A_encode=GlobalMatrixSG(A,B,HASH,Lev,Deg)
% Global Matrix construction from Matrices A and B
% by Hash Table
% This code is to use the new approach for generating global matrix
% only works for
% compression = 3
% By Lin, 11/14/2018

global hash_format

Dim = 2;


ComLev = perm_leq( Dim, Lev);
ComLevIndex = [];

for i = 1:size(ComLev,1)
    
    Lev_loc = ComLev(i,:);
    
    for d = 1:Dim
        Cell_loc{d}=[0:2^max(Lev_loc(d)-1,0)-1];
        index = LevCell2index(Lev_loc(d),Cell_loc{d});
        ComLevIndex{i}.lid{d} = index;
    end
    
    [Cel1,Cel2] = meshgrid(Cell_loc{:}); % need more work to generalize d
    
%     nz1=size(Cell_loc{1},2);
%     nz2=size(Cell_loc{2},2);
%     Cel1 = zeros(nz1*nz2,1);
%     Cel2 = zeros(nz1*nz2,1);
%     for k = 1:nz1
%     Cel1(k:nz2:end)=Cell_loc{1}(k);
%     end
%     for k = 1:nz2
%     Cel2((k-1)*nz1+1:k*nz1) = Cell_loc{2}(k);
%     end
    
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

for i = 1:size(ComLev,1)
    for j = 1:size(ComLev,1)
        IndexI = ComLevIndex{i}.gid;
        sizeI = size(IndexI,2);
        IndexJ = ComLevIndex{j}.gid;
        sizeJ = size(IndexJ,2);
        
        for d = 1:Dim
            index_I{d} = ComLevIndex{i}.lid{d};
            index_J{d} = ComLevIndex{j}.lid{d};
        end
        

%         if i == 12 || j == 12
%             111111
%         end
        
        tmpA=A( Deg*(index_I{1}(1)-1)+1:Deg*(index_I{1}(end)),...
            Deg*(index_J{1}(1)-1)+1:Deg*(index_J{1}(end)) );
        An1 = size(tmpA,1)/Deg;
        An2 = size(tmpA,2)/Deg;
        tmpB=B( Deg*(index_I{2}(1)-1)+1:Deg*(index_I{2}(end)),...
            Deg*(index_J{2}(1)-1)+1:Deg*(index_J{2}(end)) );
        Bn1 = size(tmpB,1)/Deg;
        Bn2 = size(tmpB,2)/Deg;

%         if norm(tmpA)>1e-10 && norm(tmpB)>1e-10
            
            A_encode{count}.A1=tmpA;
            A_encode{count}.A2=tmpB;
            
            IndexI = Deg^2*(IndexI(1)-1)+1:Deg^2*IndexI(end);
            IndexI = reshape(IndexI,An1*Deg,Bn1*Deg);
            t = mat2cell(IndexI,repmat(Deg,An1,1),repmat(Deg,Bn1,1));
            
            tmp_IndexI = [];
            
            for k = 1:Deg
                t = IndexI(:,k:Deg:end);
                tmp_IndexI = [tmp_IndexI;t(:)];
            end
%             t1 = IndexI(:,1:Deg:end);
%             t2 = IndexI(:,2:Deg:end);
%             IndexI = [t1(:);t2(:)];
            IndexI = tmp_IndexI;
            
            IndexJ = Deg^2*(IndexJ(1)-1)+1:Deg^2*IndexJ(end);
            IndexJ = reshape(IndexJ,An2*Deg,Bn2*Deg);
%             IndexJ = [IndexJ(:,1:Deg:end);IndexJ(:,2:Deg:end);];
%             t1 = IndexJ(:,1:Deg:end);
%             t2 = IndexJ(:,2:Deg:end);
%             IndexJ = [t1(:);t2(:)];
           tmp_IndexJ = [];
            
            for k = 1:Deg
                t = IndexJ(:,k:Deg:end);
                tmp_IndexJ = [tmp_IndexJ;t(:)];
            end
            IndexJ = tmp_IndexJ;

% %         tmpB=A( Deg*(index_I{1}(1)-1)+1:Deg*(index_I{1}(end)),...
% %             Deg*(index_J{1}(1)-1)+1:Deg*(index_J{1}(end)) );
% %         An1 = size(tmpB,1)/Deg;
% %         An2 = size(tmpB,2)/Deg;
% %         tmpA=B( Deg*(index_I{2}(1)-1)+1:Deg*(index_I{2}(end)),...
% %                 Deg*(index_J{2}(1)-1)+1:Deg*(index_J{2}(end)) );
% %         Bn1 = size(tmpA,1)/Deg;
% %         Bn2 = size(tmpA,2)/Deg;
% % 
% % 
% %         
% % %         if norm(tmpA)>1e-10 && norm(tmpB)>1e-10
% %             
% %             A_encode{count}.A1=tmpA;
% %             A_encode{count}.A2=tmpB;
% %             
% %             IndexI = Deg^2*(IndexI(1)-1)+1:Deg^2*IndexI(end);
% %             IndexI = reshape(IndexI,An1*Deg,Bn1*Deg);
% % %             t = mat2cell(IndexI,repmat(Deg,An2,1),repmat(Deg,Bn2,1));
% %             
% %             tmp_IndexI = [];
% %             
% %             for k = 1:Deg
% %                 t = IndexI(:,k:Deg:end);
% %                 tmp_IndexI = [tmp_IndexI;t(:)];
% %             end
% % 
% %             IndexI = tmp_IndexI;
% %             
% %             IndexJ = Deg^2*(IndexJ(1)-1)+1:Deg^2*IndexJ(end);
% %             IndexJ = reshape(IndexJ,An2*Deg,Bn2*Deg);
% % 
% %            tmp_IndexJ = [];
% %             
% %             for k = 1:Deg
% %                 t = IndexJ(:,k:Deg:end);
% %                 tmp_IndexJ = [tmp_IndexJ;t(:)];
% %             end
% %             IndexJ = tmp_IndexJ;
            
            A_encode{count}.IndexI = IndexI;
            A_encode{count}.IndexJ = IndexJ;
            IndexI
            IndexJ
            
            count = count+1;
%         end
    end
end

