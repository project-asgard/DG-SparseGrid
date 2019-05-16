function [elements,elementsIDX] = elementTable(pde,opts,lev,gridType)

nDims = numel(pde.dimensions);
is_sparse_grid = strcmp( opts.gridType, 'SG');

%%
% Setup element table as a collection of 2*nDims many sparse vectors to
% store the lev and cell info for each dim. 

N_max = double((uint64(2)^lev)^nDims); % This number is HUGE
% for d=1:nDims
%     elements.coords{d}.lev  = sparse(N_max,1);
%     elements.coords{d}.cell = sparse(N_max,1);
% end

elements.lev = sparse(N_max,nDims);
elements.pos = sparse(N_max,nDims);
elements.node_type = sparse(N_max,1);
% elements.idx = sparse(N_max,1);


%%
% Get the combinations of levels

if (is_sparse_grid)
   ptable = perm_leq( nDims, lev );
else
   ptable = perm_max( nDims, lev );
end

%%
% compute the number of cells for each combination of levels

ncase = size(ptable,1);
isize = zeros(1,ncase);

levels = zeros(1,nDims);
ipow = zeros(1,nDims);

for icase=1:ncase
   levels(1:nDims) = ptable(icase,1:nDims);
   ipow(1:nDims) = 2.^max(0,levels(1:nDims)-1);
   isize(icase) = prod( max(1,ipow) );
end

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

for icase=1:ncase
  istart = istartv(icase);
  iend   = istartv(icase+1)-1;

  levels(1:nDims) = ptable(icase,1:nDims);
  index_set = lev_cell_to_singleD_index_set( levels(1:nDims) );
  
  for i=istart:iend
      
     icells = index_set(i-istart+1,:);
     
     %%
     % Store the index into the element table for this element
     
     element_idx = lev_cell_to_element_index(pde,levels,icells);
     elementsIDX(i) = element_idx;
     
     %%
     % Set the lev and cell coordinate for each dim
     
     elements.lev(element_idx,:) = levels+1; % NOTE : have to start lev  index from 1 for sparse storage
     elements.pos(element_idx,:) = icells+1; % NOTE : have to start lev  index from 1 for sparse storage
     
%      for d=1:nDims    
%          elements.coords{d}.lev (element_idx) = levels(d)+1; % NOTE : have to start lev  index from 1 for sparse storage
%          elements.coords{d}.cell(element_idx) = icells(d)+1; % NOTE : have to start cell index from 1 for sparse storage   
%      end
     
     %%
     % Set the element type ("leaf" or "internal") by checking if ANY of
     % the lev coordinates are in the bottom level
     
     for d=1:nDims
         elements.node_type(element_idx) = 1; % 'internal'; % Internale nodes will not be checked for refinement.
         if elements.lev(element_idx,d) == lev+1
             elements.node_type(element_idx) = 2; % 'leaf'; % Leaf nodes are checked for refinement
         end
     end

  end
end


end
