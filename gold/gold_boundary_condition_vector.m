generate_data( @diffusion1, 'diffusion1', 'lev', 2, 'deg', 2);
generate_data( @diffusion1, 'diffusion1', 'lev', 4, 'deg', 4);
generate_data( @diffusion1, 'diffusion1', 'lev', 5, 'deg', 5);

generate_data( @mirror3, 'mirror3_case1', 'lev', 3, 'deg', 4);
generate_data( @mirror3, 'mirror3_case2', 'lev', 3, 'deg', 4);
generate_data( @mirror3, 'mirror3_case3', 'lev', 3, 'deg', 4);


function [ bcv ] = generate_data(pde_handle, pde_name, varargin)

data_dir = strcat("generated-inputs/boundary_conditions/");
root = get_root_folder();
[stat,msg] = mkdir ([root,'/gold/',char(data_dir)]);
out_format = strcat(data_dir,'vector_%s_l%d_d%d.dat');

opts = OPTS(varargin);  
opts.quiet = 1;
pde = pde_handle( opts );  

level = pde.dimensions{1}.lev;
degree = opts.deg;

[elements, elements_idx]    = hash_table_sparse_nD (pde.get_lev_vec, opts.max_lev, ...
                                                    opts.grid_type);
hash_table.elements         = elements;
hash_table.elements_idx     = elements_idx;


[~, pde.transform_blocks] = OperatorTwoScale_wavelet2(opts.deg,opts.max_lev);


TD = 0;
t = 0;
pde = get_coeff_mats(pde,opts,t,TD);

bcv = boundary_condition_vector(pde,opts,hash_table,t);

write_octave_like_output(sprintf(out_format, pde_name, level, degree), bcv );

end
