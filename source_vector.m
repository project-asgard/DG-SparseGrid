function fval = source_vector(pde, opts, hash_table, time)

% Returns the wavelet transformed source

nDims = numel(pde.dimensions);
nSources = numel(pde.sources)

%%
% Loop over the number of sources, each of which has nDims + time elements.

fval = 0;
for s=1:nSources
    for d=1:nDims
        fList{d} = forwardMWT(pde,d,pde.sources{s}{d},time);
    end
    fs_d{s}{nDims+1} = pde.sources{s}{nDims+1}(time);
    
    ft = pde.sources{s}{nDims+1}(time);
    fval = fval + combine_dimensions_D(pde, opts, fList,ft,hash_table,pde);
end

end