%% Generate gold data for C++ testing of matrix_plot_D component
data_dir = strcat("generated-inputs", "/", "matrix_plot_D", "/");
root = get_root_folder();
[stat,msg] = mkdir ([root,'/gold/',char(data_dir)]);

pde = check_pde(continuity2) 

% testing matrix_plot_D

num_dimensions = numel(pde.dimensions);

%% Construct the 1D multi-wavelet transform for each dimension.
for d=1:num_dimensions
    pde.dimensions{d}.FMWT = OperatorTwoScale(pde.deg,pde.dimensions{d}.lev);
end

for d=1:num_dimensions

  [Meval{d}, node{d}] = matrix_plot_D( pde, pde.dimensions{d});

end

out_format = strcat(data_dir, "matrix_plot_D");

for d=1:num_dimensions
write_octave_like_output(strcat(out_format, strcat("matrix_plot_D_", num2str(d), ".dat")), ...
			 full(Meval{d}));
end
