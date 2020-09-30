gen_realspace_transform( continuity1, 'continuity_1', 'lev', 8, 'deg', 7);
gen_realspace_transform( continuity2, 'continuity_2', 'lev', 7, 'deg', 6);
gen_realspace_transform( continuity3, 'continuity_3', 'lev', 6, 'deg', 5);
gen_realspace_transform( continuity6, 'continuity_6', 'lev', 2, 'deg', 3);

% Generate gold data for C++ testing of matrix_plot_D component
function [real_space] = gen_realspace_transform(pde, output_prefix, varargin)

data_dir = strcat("generated-inputs/", "transformations/", "matrix_plot_D/", output_prefix, "/");
root = get_root_folder();
[stat,msg] = mkdir ([root,'/gold/',char(data_dir)]);

runtime_defaults

pde = check_pde(pde, opts); 

num_dimensions = numel(pde.dimensions);

% Construct the 1D multi-wavelet transform for each dimension.
for d=1:num_dimensions
    pde.dimensions{d}.FMWT = OperatorTwoScale(opts.deg,pde.dimensions{d}.lev);
end

for d=1:num_dimensions

  [Meval{d}, node{d}] = matrix_plot_D( pde, opts, pde.dimensions{d});

end

out_format = strcat(data_dir, "matrix_plot_D");

for d=1:num_dimensions

%Write matrix to a file
write_octave_like_output(strcat(out_format, strcat("_", num2str(d - 1), ".dat")), ...
			 full(Meval{d}));
end
end
