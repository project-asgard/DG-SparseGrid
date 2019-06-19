function [forwardHash, inverseHash] = hash_table_nD (lev_vec, grid_type)
%-------------------------------------------------
% Generate  multi-dimension Hash Table s.t n1+n2<=Lev
% Input: Lev:: Level information
%        Dim:: Dimensionality
% Output: forwardHash:: HashTable
%         inverseHash:: Inverse Looking up for Hash
% Major Change:: ignoring the Deg from mesh
% Adding the 1D index into HashTable with
%   (Lev_1D,Cell_1D)->Index_1D
% so in 2D the inv = (Lev_1,Lev_2,Cell_1,Cell_2,Index_1,Index_2)
%              key = [Lev_1,Lev_2,Cell_1,Cell_2]
%-------------------------------------------------
idebug = 0;

if ~exist('grid_type','var') || isempty(grid_type)
    grid_type = 'SG'; % Set default gridType to SG
end
is_sparse_grid = strcmp( grid_type, 'SG');

num_dimensions = numel(lev_vec);

for d=1:num_dimensions
    assert(lev_vec(d)==lev_vec(1),'use_oldhash not supported for dimension dependent level');
end

time_perm = tic();
if (is_sparse_grid)
    ptable = perm_leq( num_dimensions, lev_vec(1) );
else
    ptable = perm_max( num_dimensions, lev_vec(1) );
end

elapsed_time_perm = toc( time_perm);
if (idebug >= 1),
    disp(sprintf('HashTable:Dim=%d,Lev=%d,gridType=%s,time %g, size=%g',...
        num_dimensions,max(lev_vec),grid_type, ...
        elapsed_time_perm,  size(ptable,1) ));
end;

global hash_format

% Specifies the number of allowable integers in the elements of the hash key
% If more are needed, i.e., > 99, then change to 'i%3.3i_'.

hash_format =  'i%04.4d';

count=1;
forwardHash = struct(); % Empty struct array
inverseHash = {}; % Empty cell array

% -------------------------------------------
% the ptable() contains the level information
% -------------------------------------------

% ----------------------------------------------------------
% compute the number of cells for each combination of levels
% ----------------------------------------------------------
ncase = size(ptable,1);
isize = zeros(1,ncase);

levels = zeros(1,num_dimensions);
ipow = zeros(1,num_dimensions);

for icase=1:ncase,
    levels(1:num_dimensions) = ptable(icase,1:num_dimensions);
    ipow(1:num_dimensions) = 2.^max(0,levels(1:num_dimensions)-1);
    isize(icase) = prod( max(1,ipow) );
end;
total_isize = sum( isize(1:ncase) );
max_isize = max( isize(1:ncase) );

if (idebug >= 1),
    disp(sprintf('HashTable:total_isize=%g', ...
        total_isize ));
end;

% ---------------------------
% prefix sum or cumulate sum
% note   matlab
%   cumsum( [2 3 4 ])
%   returns
%           [2 5 9 ]
%
%  a a  b b b  c c c c
%  istart contains
%  1    3      6       10
% ---------------------------
istartv = zeros(1,ncase+1);
istartv(1) = 1;
istartv(2:(ncase+1)) = cumsum(isize);
istartv(2:(ncase+1)) = istartv(2:(ncase+1)) + 1;

% ------------------------------------------------------------------
% Note: option whether to append the values of lev_cell_to_singleD_index(nlev,cell)
% in the inverseHash{i}
% the "key" has sufficient information to recompute this information
% ------------------------------------------------------------------
append_index_k = 1;
index_k = zeros(1,num_dimensions);

% ------------------------------
% pre-allocate temporary storage
% ------------------------------
index_set = zeros(1, max_isize );

time_hash = tic();
for icase=1:ncase,
    istart = istartv(icase);
    iend   = istartv(icase+1)-1;
        
    levels(1:num_dimensions) = ptable(icase,1:num_dimensions);
    index_set = lev_cell_to_1D_index_set( levels(1:num_dimensions) );
    for i=istart:iend,
        icells = index_set(i-istart+1,:);
        key = [levels,icells];
        
        % ----------------------------------------------------
        % TODO: check whether insertion into hash table is thread-safe
        % ----------------------------------------------------
        forwardHash.(sprintf(hash_format,key)) = i;
               
        if (append_index_k),
            for kdim=1:num_dimensions,
                index_k(kdim) = lev_cell_to_1D_index( levels(kdim), icells(kdim));
            end;
            inverseHash{i} = [key,index_k];
        else
            inverseHash{i} = [key];
        end;
        
    end
end
elapsed_hash = toc( time_hash );
if (idebug >= 1),
    disp(sprintf('HashTable: elapsed time for hashing is %g sec', ...
        elapsed_hash ));
end

end
