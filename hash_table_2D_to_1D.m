function [hash_table_1D] = hash_table_2D_to_1D(hash_table,opts)
%Converts a 2D hash table into a 1D full-grid hash_table with standard
%kronecker indexing

% TODO -- Actual 2D to 1D hash table conversion

%assert(size(hash_table.elements.lev_p1,2) == 2);

%num_basis = numel(hash_table.elements_idx);
%x_dim = 1;

%Declare now for ordering's sake
%hash_table_1D.elements.lev_p1 = [];
%hash_table_1D.elements.pos_p1 = [];

%Gather element data which gives the index of an element in each dimension 
%with respect to that dimension's full grid indexing.  
%  ( This is the same as A_data but
%    is put here for sanity's sake and reference later )
% ele_data = full(2.^(hash_table.elements.lev_p1(hash_table.elements_idx,:)) + + ...
%                 hash_table.elements.pos_p1(hash_table.elements_idx,:) - 1);
%             
% %Create hash_table_1D element_idx based on the maximum element index in the
% %x direction of the 2D hash_table.
% hash_table_1D.elements_idx = nan(1,max(ele_data(:,x_dim)));
% 
% %Get maxiumum level for sparse matrix.
% max_lev = opts.max_lev;
% 
% %Create sparse containers
% hash_table_1D.elements.lev_p1 = sparse([],[],[],2^max_lev,1,numel(hash_table_1D.elements_idx));
% hash_table_1D.elements.pos_p1 = sparse([],[],[],2^max_lev,1,numel(hash_table_1D.elements_idx));
% 
% %Populate
% count = 1;
% for i=1:num_basis
%     idx = hash_table.elements_idx(i);
%     if idx <= 2^max_lev
%         hash_table_1D.elements_idx(count) = idx;
%         hash_table_1D.elements.lev_p1(count,1) = hash_table.elements.lev_p1(idx,1);
%         hash_table_1D.elements.pos_p1(count,1) = hash_table.elements.pos_p1(idx,1);
%         count = count + 1;
%     end
% end
% hash_table_1D.elements_idx(count:end) = [];

%Get maximum level for x
max_lev = max(hash_table.elements.lev_p1(:,1))-1;

%Create full grid standard hash_table with max_lev
hash_table_1D.elements_idx = 1:2^(max_lev);
hash_table_1D.elements.lev_p1 = sparse(ceil(log2(1:2^max_lev)))'+1;
hash_table_1D.elements.pos_p1 = (1:2^max_lev)' - ...
    floor(2.^(hash_table_1D.elements.lev_p1-2));  %Time for some magic...
           

% hash_table_1D.elements.lev_p1 = hash_table.elements.lev_p1(1:2^max_lev,1);
% hash_table_1D.elements.pos_p1 = hash_table.elements.pos_p1(1:2^max_lev,1);
% hash_table_1D.elements_idx = 1:nnz(hash_table_1D.elements.pos_p1);


end

