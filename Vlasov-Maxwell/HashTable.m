function [forwardHash,inverseHash]=HashTable(maxLev,Dim)
%---------------------------------------------------------
% Algorithm 2. Creating HashTable
% Matlab Version of
% Generate Dim-dimension Hash Table 
%       s.t sum(n(1:Dim))<=maxLev
% Input: maxLev:: Level information 
%        Dim:: Dimensionality
% Output: forwardHash:: HashTable
%         inverseHash:: Inverse Looking up for Hash
%--------------------------------------------------------

forwardHash=struct();
inverseHash={};

% All combinations from choosing Dim numbers from vector [0:maxLev]
combs = permn(0:maxLev, Dim);
nLev = zeros(1,Dim);

key = zeros(1,2*Dim);
count = 1;
for i = 1:size(combs,1)
    if sum(combs(i,:))<=maxLev
        nLev = combs(i,:);
        nCell = AllCell(nLev);
        nz = size(nCell,1);
        for ii = 1:nz
            key(1:Dim) = nLev;
            key(Dim+1:end) = nCell(ii,:);
            forwardHash.(sprintf('i%g_',key))=count;
            
            index_dim = LevCell2index(nLev,nCell(ii,:));
            inverseHash{count} = [key,index_dim];

                
            count=count+1;
        end
        
        
    end
end

dof_sparse=count-1;

forwardHash.Dim=Dim;
forwardHash.dof=dof_sparse;
forwardHash.Lev=maxLev;

end





function index = LevCell2index(Lev,Cell)
%=============================================================
% for given Lev and Cell, determine Index
% Input::  Lev is 1xDim; Cell is 1xDim
% Output:: index is 1xDim
%=============================================================
dim = length(Lev);
index = zeros(1,dim);
for ii = 1:dim
   if Lev(ii) == 0
       index(ii) = 1;
   else
       index(ii)=2.^(Lev(ii)-1)+Cell(ii)+1;
   end

end

end


