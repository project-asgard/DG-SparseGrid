function [database,inv]=HashTable(Lev)
%-------------------------------------------------
% Matlab Version of
% Generate 2D Hash Table s.t n1+n2<=Lev
% Input: Lev:: Level information
%        Dim:: Dimensionality
% Output: database:: HashTable
%         inv:: Inverse Looking up for Hash
%         inv=[n1,n2,n3,p1,p2,p3,index_1,index_2,index_3]
%         Here index_i means the 1D index for (n_i,i_i)
% Major Change:: ignoring the Deg from mesh
% Added the 1D index for (Lev1D,Cell1D)
%-------------------------------------------------


count=1;
database=struct();


for n1=0:Lev
    for p1=0:max(0,2^max(0,n1-1)-1)
        
                        key=[n1,p1];
                        database.(sprintf('i%g_',key))=count;
                        
                        index_1 = LevCell2index(n1,p1);
                        
                        inv{count}=[key,index_1];
                        
%                         key=[key,index_1,index_2,index_3];
%                         database.(sprintf('i%g_',key))=count;
                        
                        count=count+1;
                   
          
        
    end
end

dof_sparse=count-1;

database.Lev=Lev;
% database.Deg=Deg;
database.dof=dof_sparse;
end


