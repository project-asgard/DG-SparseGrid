function fval = source_vector(pde, opts, hash_table, time)

% Returns the wavelet transformed source

num_dimensions  = numel(pde.dimensions);
num_sources     = numel(pde.sources);

%%
% Loop over the number of sources, each of which has nDims + time elements.

fval = 0;
for s=1:num_sources
    for d=1:num_dimensions
        fList{d} = forward_wavelet_transform(pde.deg,pde.dimensions{d}.lev,...
            pde.dimensions{d}.domainMin,pde.dimensions{d}.domainMax,...
            pde.sources{s}{d},pde.params,time);
    end
    
    time_multiplier = pde.sources{s}{num_dimensions+1}(time);
    fval = fval + combine_dimensions_D(pde.deg, fList, time_multiplier, hash_table, opts.use_oldhash);
end

end