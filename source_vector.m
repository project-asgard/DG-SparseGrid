function fval = source_vector(pde, opts, hash_table, time)

% Returns the wavelet transformed source

num_dimensions  = numel(pde.dimensions);
num_sources     = numel(pde.sources);

%%
% Loop over the number of sources, each of which has nDims + time elements.

fval = 0;
for s=1:num_sources
    for d=1:num_dimensions
        fList{d} = forward_wavelet_transform(pde,d,pde.sources{s}{d},time);
    end
    
    time_multiplier = pde.sources{s}{num_dimensions+1}(time);
    fval = fval + combine_dimensions_D(pde, opts, fList, time_multiplier, hash_table);
end

end