function [coordinates] = get_sparse_grid_coordinates(pde)
%%
% coordinates is a (num_pts x num_dims) array of real space coordinates of
% the sparse-grid element locations.

num_elems = numel(pde.elementsIDX);
num_dims  = numel(pde.dimensions);

%%
% Get center coordinates for each element

coordinates = zeros(num_elems,num_dims);
for elem=1:num_elems
    idx = pde.elementsIDX(elem);
    coordinates(elem,:) = getMyRealSpaceCoord(pde,idx);
end

%% 
% Get deg many points in each dimension for each element

% TODO


end