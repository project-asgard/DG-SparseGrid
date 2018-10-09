function [forwardHash,inverseHash,index1D] = HashTable(Lev,Dim)
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

global hash_format

% Specifies the number of allowable integers in the elements of the hash key
% If more are needed, i.e., > 99, then change to 'i%3.3i_'.

hash_format =  'i%04.4d_';

count=1;
forwardHash = struct(); % Empty struct array
inverseHash = {}; % Empty cell array

combs = perm_leq(Dim,Lev);
nLev = zeros(1,Dim);
key = zeros(1,2*Dim);
count = 1;

index1D = [];

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
        

        key = [key(1,1) key(1,2) key(2,1) key(2,2)];
        coord = [coord(1,1) coord(1,2) coord(2,1) coord(2,2) coord(3,1) coord(3,2)];
        
        keystr = getHashKeyStr(key);
        inverseHash{count} = coord;
        forwardHash.(keystr)=count;
        
                if key(1,2) == 0
            index1D = [index1D;count];
                end
        
        count=count+1;
        % for n1=0:Lev
        %     for i1=0:max(0,2^max(0,n1-1)-1)
        %
        %         for n2=0:Lev-n1
        %             for i2=0:max(0,2^max(0,n2-1)-1)
        
        %                 key=[n1,n2,i1,i2];
        %                 forwardHash.(sprintf(hash_format,key)) = count;
        %
        %                 % Linearize the heirarchial multi-index for each dimension.
        %
        %                 index_1 = LevCell2index(n1,i1);
        %                 index_2 = LevCell2index(n2,i2);
        %
        %                 inverseHash{count} = [key,index_1,index_2];
        %
        %                 count = count+1;
        %             end
        %         end
        %
        %     end
        % end
    end
end

% Add some other useful information to the forwardHash struct

forwardHash.Lev = Lev;
forwardHash.Dim = Dim;

end
