function [database,Inv]=HashTable(Lev_1,Lev_2,Deg,Dim)
%-------------------------------------------------
% Matlab Version of
% Generate 2D Hash Table s.t n1+n2<=Lev
% Input: Lev_1:: Level information for 1-component
%        Lev_2:: Level information for 2-component
%        Deg:: Degree of polynomial
%        Dim:: Dimensionality
% Output: database:: HashTable
%         Inv:: Inverse Looking up for Hash
% Update the Loop for k1
%-------------------------------------------------

count=1;
database=struct();
Inv=struct();

for n1=0:Lev_1
    for i1=0:max(0,2^max(0,n1-1)-1)
        for k1=0:Deg-1
            
            for n2=0:Lev_2-n1
                for i2=0:max(0,2^max(0,n2-1)-1)
                    for k2=0:Deg-1
                        
                        key=[n1,n2,i1,i2,k1,k2];
                        database.(sprintf('i%g_',key))=count;
                        inv{count}=key;
                        
                        % local position for v
                        if n1==0
                            index_1=Deg*i1+k1+1;
                        else
                            index_1=Deg*2^(n1-1)+Deg*i1+k1+1;
                        end
                        % local position for x
                        if n2==0
                            index_2=Deg*i2+k2+1;
                        else
                            index_2=Deg*2^(n2-1)+Deg*i2+k2+1;
                        end
                        
                        Inv.x1(count)=index_1;
                        Inv.x2(count)=index_2;
                        
                        
                        count=count+1;
                        
                        
                    end
                end
                
            end
        end
        
    end
end

dof_sparse=count-1;

database.Dim=Dim;
database.Deg=Deg;
database.dof=dof_sparse;
database.Lev_1=Lev_1;
database.Lev_2=Lev_2;
database.Lev=min(Lev_1,Lev_2);




