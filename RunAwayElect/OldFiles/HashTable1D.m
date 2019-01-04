% function [Hash,IHash,LeafHash,ILeaf,FineIndex] = HashTable1D(Lev)
function [Hash,IHash,FineIndex] = HashTable1D(Lev)
%-------------------------------------------------
% Generate 1D Hash Table
% Input: Lev:: Level information
%
% Output: forwardHash:: HashTable
%         inverseHash:: Inverse Looking up for Hash
% Major Change:: ignoring the Deg from mesh
% Adding the 1D index into HashTable with
%   (Lev_1D,Cell_1D)->Index_1D
% so the inv = (Lev,Cel,Deg,Index)
%        key = [Lev,Cel,Deg]
%-------------------------------------------------




count = 1;
count_Leaf = 1;
Hash = struct(); % Empty struct array
IHash = {}; % Empty cell array

for n1=0:Lev
    for i1=0:max(0,2^max(0,n1-1)-1)
        %         for k1=0:Deg
        
        
        
        key=[n1,i1];
        Hash.(sprintf('i%g_',key)) = count;
        
        % Linearize the heirarchial multi-index for each dimension.
        if n1==0
            index_1=i1+1;
        else
            index_1=2^(n1-1)+i1+1;
        end
        
        
        
        IHash{count} = [key,index_1];
        
        
        
       if n1 == Lev
           LeafHash.(sprintf('i%g_',key)) = count;
           ILeaf{count_Leaf} = [key,index_1];
           FineIndex(count_Leaf) = count;
           count_Leaf = count_Leaf+1;
       end
       
       count = count+1;
    end
end

% Add some other useful information to the forwardHash struct

% forwardHash.Lev = Lev;


end
