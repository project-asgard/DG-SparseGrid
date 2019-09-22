function [real_space] = gold_wavelet_to_realspace(pde, output_prefix, varargin)
% Generate gold data for C++ testing of wavelet_to_realspace component

data_dir = strcat("generated-inputs", "/", "wavelet_to_realspace", "/", output_prefix, "/");
root = get_root_folder();
[stat,msg] = mkdir ([root,'/gold/',char(data_dir)]);

runtime_defaults

opts.use_oldhash = true;

pde = check_pde(pde); 

num_dimensions = numel(pde.dimensions);

% Construct the 1D multi-wavelet transform for each dimension.
for d=1:num_dimensions
    pde.dimensions{d}.FMWT = OperatorTwoScale(pde.deg,pde.dimensions{d}.lev);
end

for d=1:num_dimensions

  [Meval{d}, node{d}] = matrix_plot_D( pde, pde.dimensions{d});

end

% create the hash table - taken directly from asgard.m
[HASH,hash_table] = hash_table_nD(pde.lev_vec, opts.grid_type);

% simple function to transform from wavelet to realspace 
domain = [ 0:((pde.deg^num_dimensions * length(hash_table) )-1) ];

%fval = 1000 * ( sin( pi * domain ) );
fval = 2 * domain;


real_space = wavelet_to_realspace( pde, opts, Meval, fval, hash_table );

% write the real space solution vector to a file
out_format = strcat(data_dir, "wavelet_to_realspace");

write_octave_like_output(strcat(out_format, ".dat"), real_space );

end
