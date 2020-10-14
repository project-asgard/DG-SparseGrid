function [f_realspace_nD_fixed_grid_nodups,nodes_nodups]...
    = get_realspace_fixed_grid_solution(pde,opts,fval,hash_table)

num_fixed_grid = 256;

num_dimensions = numel(pde.dimensions);

if nargin<4
    if opts.use_oldhash
        [HASH,hash_table] = hash_table_nD(pde.lev_vec, opts.grid_type);
    else
        [elements, elements_idx]    = hash_table_sparse_nD (pde.lev_vec, pde.max_lev, opts.grid_type);
        hash_table.elements         = elements;
        hash_table.elements_idx     = elements_idx; % only to get the same order as the old hash table
    end
end

for d=1:num_dimensions
    nodes_nodups{d} = linspace(pde.dimensions{d}.min,pde.dimensions{d}.max,num_fixed_grid);
    [Meval{d},nodes{d},nodes_count{d}] = matrix_plot_D(pde,opts,pde.dimensions{d},nodes_nodups{d});
end

fval_realspace_fixed_grid = wavelet_to_realspace(pde,opts,Meval,fval,hash_table);

f_realspace_nD_fixed_grid = singleD_to_multiD(num_dimensions,fval_realspace_fixed_grid,nodes);
f_realspace_nD_fixed_grid_nodups = ...
    remove_duplicates(num_dimensions,f_realspace_nD_fixed_grid,nodes_nodups,nodes_count);

end
