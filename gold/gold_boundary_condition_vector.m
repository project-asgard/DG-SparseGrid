generate_data( diffusion1, 'diffusion1', 'lev', 2, 'deg', 2);
generate_data( diffusion1, 'diffusion1', 'lev', 4, 'deg', 4);
generate_data( diffusion1, 'diffusion1', 'lev', 5, 'deg', 5);

function [ bcv ] = generate_data(pde, output_prefix, varargin)

data_dir = strcat("generated-inputs", "/", "boundary_condition_vector", "/", output_prefix, "/");
root = get_root_folder();
[stat,msg] = mkdir ([root,'/gold/',char(data_dir)]);

runtime_defaults

opts.use_oldhash = true;

pde = check_pde(pde,opts);

level = pde.dimensions{1}.lev;
degree = pde.deg;

% create the hash table - taken directly from asgard.m
[HASH,hash_table] = hash_table_nD(pde.lev_vec, opts.grid_type);

for d=1:num_dimensions
    pde.dimensions{d}.FMWT = OperatorTwoScale(pde.deg,pde.dimensions{d}.lev);
end


t = 0;
TD = 0;
pde = get_coeff_mats(pde,t,TD);

TD = 1;
pde = get_coeff_mats(pde,t,TD);
        
bcv = boundary_condition_vector(pde,opts,hash_table,t);

out_format = strcat(data_dir,...
                    "boundary_condition_vector",...
                    "_l", num2str( level ),...
                    "_d", num2str( degree ) );

write_octave_like_output(strcat(out_format, ".dat"), bcv );

end
