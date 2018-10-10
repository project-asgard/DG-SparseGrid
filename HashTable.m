function [forwardHash,inverseHash,hashIndex1,hashIndex2] = HashTable(Lev,Dim)
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
%-------------------------------------------------

count1=1;
forwardHash = struct(); % Empty struct array
inverseHash = {}; % Empty cell array

combs = perm_leq(Dim,Lev);
nLev = zeros(1,Dim);
key = zeros(1,2*Dim);

count1 = 1;
count2 = 1;
hashIndex1 = [];
hashIndex2 = [];

for i = 1:size(combs,1)
    nLev = combs(i,:);
    for d = 1:Dim
        nMax = max(0,2^max(0,nLev(d)-1)-1);
        if nMax ==0
            value{d} = [0];
        else
            value{d} = [0:nMax];
        end
    end
    nCell = allcomb(value{:});
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
        
%         % Remap the new hash layout to the old one
%         key = [key(1,1) key(1,2) key(2,1) key(2,2)];
%         coord = [coord(1,1) coord(1,2) coord(2,1) coord(2,2) coord(3,1) coord(3,2)];
        
        keystr = getHashKeyStr(key);
        inverseHash{count1} = coord;
        forwardHash.(keystr)=count1;
        
        % Construct hash indices for each dimension
        if key(1,1) == 0
            hashIndex1 = [hashIndex1;count1];
        end
        if key(1,2) == 0
            hashIndex2 = [hashIndex2;count2];
        end
        
        count1=count1+1;
        count2=count2+1;
        
    end
end


% Add some other useful information to the forwardHash struct

forwardHash.Lev = Lev;
forwardHash.Dim = Dim;

end
