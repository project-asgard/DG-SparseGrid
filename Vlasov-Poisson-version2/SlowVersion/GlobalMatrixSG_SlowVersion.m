function A_encode=GlobalMatrixSG_SlowVersion(A,B,HASHInv,Con2D,Deg)
% Global Matrix construction from Matrices A and B
% by Hash Table
% Expanding Deg
% w.r.t HASH table
%-----------------------------------------------------------
HASHDOF = size(HASHInv,2);

count = 1;

for i = 1:HASHDOF
    % get the mesh information for Row i with 
    % (Lev1,Lev2,Cel1,Cel2,Index1,Index2)
    ll = HASHInv{i};
    
    % 1D index through (Lev1,Cel1) and (Lev2,Cel2)
    I1 = ll(5);
    I2 = ll(6);
    
    % the nonzero columns for Row i
    j = Con2D{i};
    
    % the key for nonzero columns
    JJ = [HASHInv{j}];
    
    % 1D indices for (Lev_1D,Cell_1D)
    J1 = JJ(5:6:end);
    J2 = JJ(6:6:end);
    
    % loop of Deg
    for k1 = 1:Deg
        index_I1 = (I1-1)*Deg+k1;
        for k2 = 1:Deg
            
            index_I2 = (I2-1)*Deg+k2;
            IndexI = Deg^2*(i-1)+Deg*(k1-1)+k2;
            
            for jjj = 1:size(J1,2)
                
                for kk1 = 1:Deg
                    
                    index_J1 = (J1(jjj)-1)*Deg+kk1;
                    
                    for kk2 = 1:Deg
                        
                        
                        index_J2 = (J2(jjj)-1)*Deg+kk2;
                        IndexJ=Deg^2*(j(jjj)-1)+Deg*(kk1-1)+kk2;
                        
                        tmpA=A(index_I1,index_J1);
                        tmpB=B(index_I2,index_J2);
                        
                        A_encode{count}.IndexI = IndexI;
                        
                        A_encode{count}.A1=tmpA;
                        A_encode{count}.A2=tmpB;
                        
                        
                        A_encode{count}.IndexJ = IndexJ;
                        
                        count = count+1;
                    end
                end
                
            end
            
        end
    end
    
    
    
    
    
    
    
end
