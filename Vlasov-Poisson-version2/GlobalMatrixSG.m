function A_encode=GlobalMatrixSG(A,B,HASHInv,Con2D,Deg)
% Global Matrix construction from Matrices A and B
% by Hash Table
% 
% w.r.t HASH table
%-----------------------------------------------------------
HASHDOF = size(HASHInv,2);

count=1;

for i=1:HASHDOF
    % get the mesh information for Row i with 
    % (Lev1,Lev2,Cel1,Cel2,Index1,Index2)
    ll=HASHInv{i};

    % 1D index through (Lev1,Cel1) and (Lev2,Cel2)
    I1 = ll(5);
    I2 = ll(6);
    
    % the nonzero columns for Row i
    j=Con2D{i};
    
    % the key for nonzero columns
    JJ=[HASHInv{j}];
    
    % 1D indices for (Lev_1D,Cell_1D)
    J1=JJ(5:6:end);
    J2=JJ(6:6:end);
    

    % A(index_I1,index_J1) and B(index_I2,index_J2)
    %----------------------
    % More work can be done but not work now
    %----------------------
%     [~,JJ1,~]=unique(J1);

   
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


