function [coordinates] = get_sparse_grid_coordinates(pde, opts, hash_table)
%%
% coordinates is a (num_pts x num_dims) array of real space coordinates of
% the sparse-grid element locations.

if opts.use_oldhash
    num_elems = numel(hash_table);  
else
    num_elems = numel(hash_table.elements_idx);
end
num_dims  = numel(pde.dimensions);

%%
% Get center coordinates for each element

coordinates = zeros(num_elems,num_dims);
for elem=1:num_elems
    idx = hash_table.elements_idx(elem);
    coordinates(elem,:) = getMyRealSpaceCoord(pde, opts, hash_table, idx);
end

%% 
% Get deg many points in each dimension for each element

% TODO


end