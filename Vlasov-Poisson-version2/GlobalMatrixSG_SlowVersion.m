function A_encode=GlobalMatrixSG_SlowVersion(A,B,HASH,HASHInv,Con2D,Deg)
% Global Matrix construction from Matrices A and B
% by Hash Table
% Expanding Deg
% w.r.t HASH table
%-----------------------------------------------------------
count=1;

for i=1:HASH.dof
    % get the key for Row i with (Lev1,Lev2,Cel1,Cel2)
    ll=HASHInv{i};
    n1=ll(1);c1=ll(3);
    n2=ll(2);c2=ll(4);
    
    % 1D index through (Lev1,Cel1)
    % 1D index through (Lev2,Cel2)
    I1=LevCell2index(n1,c1);
    I2=LevCell2index(n2,c2);
    
    % the nonzero columns for Row i
    j=Con2D{i};
    
    % the key for nonzero columns
    JJ=[HASHInv{j}];
    
    % (Lev,Cell) informations for component 1 and 2
    % (m1,p1) and (m2,p2)
    m1=JJ(1:4:end);p1=JJ(3:4:end);
    m2=JJ(2:4:end);p2=JJ(4:4:end);
    
    %================================================================
    % 1D index
    % This is to compute the 1D index from (Lev_1D,Cell_1D)-->index
    % index={ 2^(Lev-1)+Cell+1, if Lev~= 0
    %       { 1               , if Lev== 0
    %================================================================
    J1=LevCell2index(m1,p1);
    J2=LevCell2index(m2,p2);
    
    % loop of Deg
        for k2 = 1:Deg
            index_I2=(I2-1)*Deg+k2;
            for k1 = 1:Deg
                
                index_I1=(I1-1)*Deg+k1;
                IndexI = Deg^2*(i-1)+Deg*(k2-1)+k1;
                
                for jjj = 1:size(J1,2)
                    
                    for kk1=1:Deg
                        for kk2=1:Deg
                        
                            index_J1=(J1(jjj)-1)*Deg+kk1;
                            index_J2=(J2(jjj)-1)*Deg+kk2;

                            tmpA=A(index_I1,index_J1);
                            tmpB=B(index_I2,index_J2);
        
                            A_encode{count}.IndexI = IndexI;
    
                            A_encode{count}.A1=tmpA;
                            A_encode{count}.A2=tmpB;

                            IndexJ=Deg^2*(j(jjj)-1)+Deg*(kk2-1)+kk1;
                            A_encode{count}.IndexJ = IndexJ;
          
                            count = count+1;
                        end
                    end
        
                end 
                
            end
        end
            
    
    
   
    
    

end
