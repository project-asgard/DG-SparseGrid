function fval = md_eval_function(opts, dims, params, md_funcs, hash_table, t)

% Returns a multi-D wavelet space function defined by "md_funcs"

num_dimensions  = numel(dims);
num_funcs       = numel(md_funcs);

%%
% Loop over the number of sources, each of which has nDims + time elements.

fval = 0;
for s=1:num_funcs
    for d=1:num_dimensions
        fList{d} = forward_wavelet_transform(deg,dims{d}.lev, ...
            dims{d}.domainMin,dims{d}.domainMax, ...
            md_funcs{s}{d},params,t);
    end
    
    time_multiplier = md_funcs{s}{num_dimensions+1}(t);
    fval = fval + combine_dimensions_D(deg, fList, time_multiplier, hash_table, opts.use_oldhash);
end

end