generate_data( diffusion1, 'diffusion1', 'lev', 2, 'deg', 2);
generate_data( diffusion1, 'diffusion1', 'lev', 4, 'deg', 4);
generate_data( diffusion1, 'diffusion1', 'lev', 5, 'deg', 5);

function [ bcv ] = generate_data(pde, pde_name, varargin)

data_dir = strcat("generated-inputs/boundary_conditions/");
root = get_root_folder();
[stat,msg] = mkdir ([root,'/gold/',char(data_dir)]);
out_format = strcat(data_dir,'vector_%s_l%d_d%d.dat');

runtime_defaults

pde = check_pde(pde,opts);

level = pde.dimensions{1}.lev;
degree = pde.deg;

[elements, elements_idx]    = hash_table_sparse_nD (pde.lev_vec, opts.max_lev, ...
                                                    opts.grid_type);
hash_table.elements         = elements;
hash_table.elements_idx     = elements_idx;


for d=1:num_dimensions
    pde.dimensions{d}.FMWT = OperatorTwoScale(pde.deg,pde.dimensions{d}.lev);
end

TD = 0;
t = 0;
pde = get_coeff_mats(pde,t,TD);

bcv = boundary_condition_vector(pde,opts,hash_table,t);

write_octave_like_output(sprintf(out_format, pde_name, level, degree), bcv );

end