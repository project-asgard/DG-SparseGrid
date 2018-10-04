function [forwardHash,inverseHash] = HashTable2(Lev,Dim)
%-------------------------------------------------
% Generate 2D Hash Table s.t n1+n2<=Lev
% Input: Lev:: Level information
%        Dim:: Dimensionality
% Output: forwardHash:: HashTable
%         inverseHash:: Inverse Looking up for Hash
% Major Change:: ignoring the Deg from mesh
% Adding the 1D index into HashTable with
%   (Lev_1D,Cell_1D)->Index_1D
% so the inv = (Lev_1,Lev_2,Cell_1,Cell_2,Index_1,Index_2)
%        key = [Lev_1,Lev_2,Cell_1,Cell_2]
% Note:: 10/02--the structure of Hash table and Inverse Hash have been
% modified
%-------------------------------------------------

count=1;
forwardHash = struct(); % Empty struct array
inverseHash = {}; % Empty cell array

combs = perm_leq(Dim,Lev);

nLev = zeros(1,Dim);

key = zeros(1,2*Dim);
count = 1;

for i = 1:size(combs,1)
    
    nLev = combs(i,:);
    nCell = AllCell(nLev);
    nz = size(nCell,1);
    
    for ii = 1:nz
        
        key = [];
        coord = [];
        
        for d = 1:Dim
            
            lev = nLev(d);
            cell = nCell(ii,d);
            indx = LevCell2index(lev,cell);
            
            key(1:2,d) =[lev,cell];
            
            coord(1:3,d) = [lev,cell,indx];
            
        end
        
        keystr = getHashKeyStr(key);
        inverseHash{count} = coord;
        forwardHash.(keystr)=count;
        count=count+1;
    end
end

dof_sparse=count-1;

forwardHash.Dim=Dim;
forwardHash.dof=dof_sparse;
forwardHash.Lev=Lev;

end
