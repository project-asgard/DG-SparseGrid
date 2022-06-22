function [f] = fill_from_hash_table(pde,opts,hash_table_old,hash_table,f_old)
%After refinement, the function needs to be backfilled. 
%This algorithm does that.  

%Assert that the new hash_table is just a refinement of the old one
num_ele_old = numel(hash_table_old.elements_idx);
num_ele = numel(hash_table.elements_idx);

num_dims = numel(pde.dimensions);
deg = opts.deg;

cell_dof = deg^num_dims;

assert(sum(hash_table_old.elements_idx-hash_table.elements_idx(1:num_ele_old)) == 0) %This is an integer sum so setting to zero is okay.
assert(num_ele_old*cell_dof == numel(f_old));

f = zeros(num_ele*cell_dof,1);
f(1:num_ele_old*cell_dof) = f_old;


end

