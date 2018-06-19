function [database,Inv]=HashTable(Lev,Dim)
%-------------------------------------------------
% Matlab Version of
% Generate 2D Hash Table s.t n1+n2<=Lev
% Input: Lev_1:: Level information for 1-component
%        Lev_2:: Level information for 2-component
%        Deg:: Degree of polynomial
%        Dim:: Dimensionality
% Output: database:: HashTable
%         Inv:: Inverse Looking up for Hash
%-------------------------------------------------

count=1;
database=struct();
Inv=struct();




%    K = 3; % The required sum
%    n = 4;  % The number of elements in the rows
%    c = nmultichoosek(values, k);
%    m = size(c,1);
%    A = zeros(m,n);
%    for ix = 1:m
%      A(ix,:) = diff([0,c(ix,:),K+1]);
%    end
   
% All combination
combs = nmultichoosek(0:Lev, Dim)

key = zeros(1,2*Dim);
count = 1;
for i = 1:size(combs,1)
    if sum(combs(i,:))<Lev
        for loc_dim = 1:Dim
            nlev = combs(i,loc_dim);
            for loc_cell = 0:max(0,2^max(0,nlev-1)-1)
                
                key(loc_dim)   = nlev;
                key(Dim+loc_dim) = loc_cell;
                key
                
                [loc_dim,loc_cell,count]
                database.(sprintf('i%g_',key))=count;
                
                count=count+1;
            end
        end
    end
end


dof_sparse=count-1;

database.Dim=Dim;
% database.Deg=Deg;
database.dof=dof_sparse;
% database.Lev_1=Lev_1;
% database.Lev_2=Lev_2;
% database.Lev=min(Lev_1,Lev_2);

function combs = nmultichoosek(values, k)
%// Return number of multisubsets or actual multisubsets.
if numel(values)==1
    n = values;
    combs = nchoosek(n+k-1,k);
else
    n = numel(values);
    combs = bsxfun(@minus, nchoosek(1:n+k-1,k), 0:k-1);
    combs = reshape(values(combs),[],k);
end






