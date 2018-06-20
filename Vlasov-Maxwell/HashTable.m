function [forwardHash,inverseHash]=HashTable(Lev,Dim)
%---------------------------------------------------------
% Matlab Version of
% Generate Dim-dimension Hash Table s.t sum(n(1:Dim))<=Lev
% Input: Lev:: Level information 
%        Dim:: Dimensionality
% Output: forwardHash:: HashTable
%         inverseHash:: Inverse Looking up for Hash
%--------------------------------------------------------

forwardHash=struct();
inverseHash={};

% All combinations from choosing Dim numbers from vector [0:Lev]
combs = permn(0:Lev, Dim);
nLev = zeros(1,Dim);

key = zeros(1,2*Dim);
count = 1;
for i = 1:size(combs,1)
    if sum(combs(i,:))<=Lev
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
forwardHash.Lev=Lev;

end

function [M, I] = permn(V, N)
%========================================================
% Compute all N combinators with repetition from vector V
% PERMN - permutations with repetition
%========================================================

nV = numel(V) ;

if nV==0 || N == 0
    M = zeros(nV,N) ;
    I = zeros(nV,N) ;
    
elseif N == 1
    % return column vectors
    M = V(:) ;
    I = (1:nV).' ;
else
    % this is faster than the math trick used for the call with three
    % arguments.
    [Y{N:-1:1}] = ndgrid(1:nV) ;
    I = reshape(cat(N+1,Y{:}),[],N) ;
    M = V(I) ;
end
end


function cell = AllCell(Lev)
%===================================================
% Compute all the cell information from Lev
% Lev is 1xDim array
% cell is *xDim array, and each row corresponding to
%       varying Lev(i), i=1~Dim
%===================================================
dim = length(Lev);

nMax = zeros(1,dim);
for ii =1:dim
    nMax(ii) = max(0,2^max(0,Lev(ii)-1)-1);
end
nz = prod(nMax+1);
    

cell = zeros(nz,dim);
for ii=1:dim
    % Both of the following two methods can work
    % Method 1
    % tmp=repmat([0:nMax(ii)]',nz/(nMax(ii)+1),1);
    % cell(:,ii)=tmp(:);
    % Method 2
    cell(:,ii)=repmat([0:nMax(ii)]',nz/(nMax(ii)+1),1);
end

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






