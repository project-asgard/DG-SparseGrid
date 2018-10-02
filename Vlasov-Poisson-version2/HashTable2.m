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
%-------------------------------------------------

global hash_format

% Specifies the number of allowable integers in the elements of the hash key
% If more are needed, i.e., > 99, then change to 'i%3.3i_'.

hash_format =  'i%04.4d_';

count=1;
forwardHash = struct(); % Empty struct array
inverseHash = {}; % Empty cell array

% result = perm_leq(Dim,Lev);
% icount = perm_leq_count(Dim,Lev);

combs = perm_leq(Dim,Lev);
%icount = perm_leq_count(Dim,maxLev);

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
%         key(1:Dim) = nLev;
%         key(Dim+1:end) = nCell(ii,:);
%         forwardHash.(sprintf(hash_format,key))=count;
%         
%         index_dim = LevCell2index(nLev,nCell(ii,:));
% %         inverseHash{count} = [key,index_dim];
%         
%         inverseHash{count} = [key,index_dim];
        lev = nLev(d);
        cell = nCell(ii,d);
        indx = LevCell2index(lev,cell);
%         key((d-1)*2+1) = lev;%nLev;
%         key(Dim+1:end) = nCell(ii,:);
        key(1:2,d) =[lev,cell];
        
%         forwardHash.(sprintf(hash_format,key))=count;
        
%         index_dim = LevCell2index(nLev,nCell(ii,:));
%         inverseHash{count} = [key,index_dim];
        coord(1:3,d) = [lev,cell,indx]; 
        
        end
        keystr = sprintf(hash_format,key(:));
        inverseHash{count} = coord;%[key,index_dim];
        forwardHash.(keystr)=count;
        count=count+1;
    end
    
end

dof_sparse=count-1;

forwardHash.Dim=Dim;
forwardHash.dof=dof_sparse;
forwardHash.Lev=Lev;


end
