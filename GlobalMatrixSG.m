function A_encode=GlobalMatrixSG(A,B,HASH,Lev,Deg)
% Global Matrix construction from Matrices A and B
% by Hash Table
% This code is to use the new approach for generating global matrix
% only works for
% compression = 3
% By Lin, 11/14/2018
% 
% w.r.t HASH table
%-----------------------------------------------------------
% % HASHDOF = size(HASHInv,2);
% % 
% % count=1;
% % 
% % for i=1:HASHDOF
% %     % get the mesh information for Row i with 
% %     % (Lev1,Lev2,Cel1,Cel2,Index1,Index2)
% %     ll=HASHInv{i};
% % 
% %     % 1D index through (Lev1,Cel1) and (Lev2,Cel2)
% %     I1 = ll(5);
% %     I2 = ll(6);
% %     
% %     % the nonzero columns for Row i
% %     j=Con2D{i};
% %     
% %     % the key for nonzero columns
% %     JJ=[HASHInv{j}];
% %     
% %     % 1D indices for (Lev_1D,Cell_1D)
% %     J1=JJ(5:6:end);
% %     J2=JJ(6:6:end);
% %     
% % 
% %     % A(index_I1,index_J1) and B(index_I2,index_J2)
% %     %----------------------
% %     % More work can be done but not work now
% %     %----------------------
% % %     [~,JJ1,~]=unique(J1);
% % 
% %    
% %     index_I1=[(I1-1)*Deg+1:I1*Deg];
% %     index_I2=[(I2-1)*Deg+1:I2*Deg];
% %     
% %     %================================================================
% %     % Since the HashTable only contains the (Lev_2D,Cell_2D), we have
% %     % to associate the 2D global index with Deg::
% %     % For example, if the HashTable with index value "1"
% %     % Then the global index (with Deg) will be [1,2,3,4,...,Deg^2]
% %     %================================================================
% %     IndexI=[Deg^2*(i-1)+1:Deg^2*i];
% %     
% %     for jjj = 1:size(J1,2)
% %         
% %         index_J1=[(J1(jjj)-1)*Deg+1:J1(jjj)*Deg];
% %         index_J2=[(J2(jjj)-1)*Deg+1:J2(jjj)*Deg];
% % 
% %         tmpA=A(index_I1,index_J1);
% %         tmpB=B(index_I2,index_J2);
% %         
% %         A_encode{count}.IndexI = IndexI;
% %     
% %         A_encode{count}.A1=tmpA;
% %         A_encode{count}.A2=tmpB;
% % 
% %         IndexJ=[Deg^2*(j(jjj)-1)+1:Deg^2*j(jjj)];
% %         A_encode{count}.IndexJ = IndexJ;
% %           
% %         count = count+1;
% %         
% %     end 
% % end

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
        IndexJ = ComLevIndex{j}.gid;
        
        for d = 1:Dim
            index_I{d} = ComLevIndex{i}.lid{d};
            index_J{d} = ComLevIndex{j}.lid{d};
        end
        
        tmpA=A( Deg*(index_I{1}(1)-1)+1:Deg*(index_I{1}(end)),...
                Deg*(index_J{1}(1)-1)+1:Deg*(index_J{1}(end)) );
        tmpB=B( Deg*(index_I{2}(1)-1)+1:Deg*(index_I{2}(end)),...
                Deg*(index_J{2}(1)-1)+1:Deg*(index_J{2}(end)) ); 
        
        
    
        A_encode{count}.A1=tmpA;
        A_encode{count}.A2=tmpB;

        IndexI = Deg^2*(IndexI(1)-1)+1:Deg^2*IndexI(end);
        IndexJ = Deg^2*(IndexJ(1)-1)+1:Deg^2*IndexJ(end);
        
        IndexI = reshape(IndexI,numel(IndexI)/Deg,Deg);
%         IndexI = reshape(IndexI,Deg,numel(IndexI)/Deg);
        IndexI = IndexI';
%         IndexJ = reshape(IndexJ,Deg,numel(IndexJ)/Deg);
        IndexJ = reshape(IndexJ,numel(IndexJ)/Deg,Deg);
        IndexJ = IndexJ';
        
        A_encode{count}.IndexI = IndexI(:);
        A_encode{count}.IndexJ = IndexJ(:);
        
        count = count+1;
    end    
end

