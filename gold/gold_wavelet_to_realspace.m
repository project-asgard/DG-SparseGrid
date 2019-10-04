function [real_space] = gold_wavelet_to_realspace(pde, varargin)
% Generate gold data for C++ testing of wavelet_to_realspace component

data_dir = strcat("generated-inputs", "/", "wavelet_to_realspace", "/");
root = get_root_folder();
[stat,msg] = mkdir ([root,'/gold/',char(data_dir)]);

runtime_defaults

pde = check_pde(continuity2); 

    %Captain! Test code 
    pde.lev_vec
    pde.max_lev
    opts.grid_type
    %Captain! End test code 

num_dimensions = numel(pde.dimensions);

% Construct the 1D multi-wavelet transform for each dimension.
for d=1:num_dimensions
    pde.dimensions{d}.FMWT = OperatorTwoScale(pde.deg,pde.dimensions{d}.lev);
end

for d=1:num_dimensions

  [Meval{d}, node{d}] = matrix_plot_D( pde, pde.dimensions{d});

end

% set only the fields needed by this test - not a full options test
%opts.use_oldhash = false;

% should this be full grid for the test? Change the values and find out
%opts.grid_type = 'FG';

% create the hash table - taken directly from asgard.m
[elements, elements_idx]    = element_table (pde.lev_vec, pde.max_lev, opts.grid_type);
hash_table.elements         = elements;
hash_table.elements_idx     = elements_idx; 

fval = initial_condition_vector(pde, opts, hash_table, 0);

fval = [ 0:(length(fval)-1) ] .* 2.0;

length(fval)

real_space = wavelet_to_realspace( pde, opts, Meval, fval, hash_table );

%length(real_space)

%real_space( 1:20 )
end

%modify fval to be something other than zero 
% Simple function in wavespace to transform into real - not a realistic example
% but hopefully useful for testing

% what pieces do I need? I need options... but just one is used right.
% done pde - from check_pde 
% done opts - only thing used is use_oldhash or nah. Set it to default and use that. also gridtype
% done meval - get that from matrix_plot_D
% fval - get that from initial_condition_vector() function
% done hash table - use element table and grid type from opts
