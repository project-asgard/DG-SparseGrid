generate_data( continuity1, 'continuity_1', 'lev', 8, 'deg', 7);
generate_data( continuity2, 'continuity_2', 'lev', 4, 'deg', 5);
generate_data( continuity3, 'continuity_3', 'lev', 3, 'deg', 4);

% Generate gold data for C++ testing of wavelet_to_realspace component
function [real_space] = generate_data(pde, pde_name, varargin)

data_dir = "generated-inputs/transformations/";
out_format = strcat(data_dir,'wavelet_to_realspace_%s.dat');
root = get_root_folder();
[stat,msg] = mkdir ([root,'/gold/',char(data_dir)]);

runtime_defaults

opts.use_oldhash = false;

pde = check_pde(pde, opts); 

num_dimensions = numel(pde.dimensions);

[~, pde.transform_blocks] = OperatorTwoScale_wavelet2(opts.deg, opts.max_lev);


for d=1:num_dimensions

  [Meval{d}, node{d}] = matrix_plot_D( pde, opts, pde.dimensions{d});

end

% create the hash table - taken directly from asgard.m
[elements, elements_idx]    = hash_table_sparse_nD (opts.lev_vec, ...
                                    opts.max_lev, opts.grid_type);
hash_table.elements         = elements;
hash_table.elements_idx     = elements_idx;

% simple function to transform from wavelet to realspace 
domain = [ 0:((opts.deg^num_dimensions * size(elements_idx,2) )-1) ];

fval = 2 * domain;


real_space = wavelet_to_realspace( pde, opts, Meval, fval, hash_table );

% write the real space solution vector to a file
write_octave_like_output(sprintf(out_format,pde_name), real_space );

end
