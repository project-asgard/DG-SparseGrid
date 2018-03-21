function [database,inv,Inv]=HashTable(Lev,Dim)
%-------------------------------------------------
% Matlab Version of
% Generate 2D Hash Table s.t n1+n2<=Lev
% Input: Lev:: Level information
%        Dim:: Dimensionality
% Output: database:: HashTable
%         inv:: Inverse Looking up for Hash
% Major Change:: ignoring the Deg from mesh
%-------------------------------------------------

global hash_format
% Specifies the number of allowable integers in the elements of the hash key
% If more are needed, i.e., > 99, then change to 'i%3.3i_'.
hash_format = 'i%2.2i_';

count=1;
database=struct();
Inv=struct();


for n1=0:Lev
    for i1=0:max(0,2^max(0,n1-1)-1)
        
        for n2=0:Lev-n1
            for i2=0:max(0,2^max(0,n2-1)-1)
                
                key=[n1,n2,i1,i2];
                database.(sprintf('i%g_',key))=count;
                
                inv{count}=key;
                
                Inv.x1(count)=LevCell2index(n1,i1);
                Inv.x2(count)=LevCell2index(n2,i2);
                
                count=count+1;
            end
        end
        
    end
end

dof_sparse=count-1;

database.Lev=Lev;
database.Dim=Dim;
database.dof=dof_sparse;
