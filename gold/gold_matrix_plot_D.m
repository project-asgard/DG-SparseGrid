% Generate gold data for C++ testing of matrix_plot_D component
function [real_space] = gold_matrix_plot_D(pde, output_prefix, varargin)

data_dir = strcat("generated-inputs", "/", "matrix_plot_D", "/", output_prefix, "/");
root = get_root_folder();
[stat,msg] = mkdir ([root,'/gold/',char(data_dir)]);

runtime_defaults

pde = check_pde(pde); 

% testing matrix_plot_D

num_dimensions = numel(pde.dimensions);

% Construct the 1D multi-wavelet transform for each dimension.
for d=1:num_dimensions
    pde.dimensions{d}.FMWT = OperatorTwoScale(pde.deg,pde.dimensions{d}.lev);
end

for d=1:num_dimensions

  [Meval{d}, node{d}] = matrix_plot_D( pde, pde.dimensions{d});

end

out_format = strcat(data_dir, "matrix_plot_D");

for d=1:num_dimensions

%Write matrix to a file
write_octave_like_output(strcat(out_format, strcat("_", num2str(d - 1), ".dat")), ...
			 full(Meval{d}));
end
