function A_encode=GlobalMatrixSG(mat,Hash)
%=========================================
% Algorithm 5: Global Assembling
%=========================================
% save the submatrix and apply A_encode to vector f
count=0;
%=====================================================
% The Index is grouped by the sum of Levels in V and X
% {I_k} together with Cell and Deg information will
%         determine the position of matrix A
% {J_k} together with Cell and Deg information will
%         determine the position of matrix B
% Note: SG needs sum_k I_k<=Lev and sum_k J_k<=Lev
% Here k = 1,2,...,Dim
%=====================================================
Deg = Hash.Deg;
Lev = Hash.Lev;
HASHDOF = size(HASHInv,2);

count=1;

for i=1:HASHDOF
    % get the mesh information for Row i with 
    % (Lev1,Lev2,Cel1,Cel2,Index1,Index2)
    ll=HASHInv{i};
    
    Ind = ones(1,6);
    % get index for each dimension: here we seperate
    % dimension x and dimension v
    Ind(1:DimX) = ll(end-Dim+1:end-DimV);
    Ind(4:3+DimV) = ll(end-DimV+1:end);
    
    for j = 1:DimX
        jind = Con1D{Ind(j)} ;
        loc_X_mat{j} = mat{j}();
    end
    for j=1:DimV
        jind = Con1D{Ind(j+3)} ;
        loc_V_mat{j} = mat{3+j}();
    end
    
    
    % find the inverse Hash from table
    
    % the nonzero columns for Row i
    j=Con2D{i};
    
    % the key for nonzero columns
    JJ=[HASHInv{j}];
    
    % 1D indices for (Lev_1D,Cell_1D)
    J1=JJ(5:6:end);
    J2=JJ(6:6:end);
    
    index_I1=[(I1-1)*Deg+1:I1*Deg];
    index_I2=[(I2-1)*Deg+1:I2*Deg];
    
    %================================================================
    % Since the HashTable only contains the (Lev_2D,Cell_2D), we have
    % to associate the 2D global index with Deg::
    % For example, if the HashTable with index value "1"
    % Then the global index (with Deg) will be [1,2,3,4,...,Deg^2]
    %================================================================
    IndexI=[Deg^2*(i-1)+1:Deg^2*i];
    
    for jjj = 1:size(J1,2)
        
        index_J1=[(J1(jjj)-1)*Deg+1:J1(jjj)*Deg];
        index_J2=[(J2(jjj)-1)*Deg+1:J2(jjj)*Deg];

        tmpA=A(index_I1,index_J1);
        tmpB=B(index_I2,index_J2);
        
        A_encode{count}.IndexI = IndexI;
    
        A_encode{count}.A1=tmpA;
        A_encode{count}.A2=tmpB;

        IndexJ=[Deg^2*(j(jjj)-1)+1:Deg^2*j(jjj)];
        A_encode{count}.IndexJ = IndexJ;
          
        count = count+1;
        
    end 
    

end


end
