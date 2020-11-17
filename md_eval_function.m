function fval = md_eval_function(opts, deg, dims, params, md_funcs, hash_table, blocks, t)

% Returns a multi-D wavelet space function defined by "md_funcs"

num_dims  = numel(dims);
num_md_funcs = numel(md_funcs);

%%
% Loop over the number of sources, each of which has nDims + time elements.

fval = 0;
for s=1:num_md_funcs
    
    md_func = md_funcs{s};
    
    for d=1:num_dims
        fList{d} = forward_wavelet_transform(deg,dims{d}.lev, ...
            dims{d}.min,dims{d}.max, ...
            md_func{d},params, blocks, t);
    end
    
    time_multiplier = md_funcs{s}{num_dims+1}(t);
    fval = fval + combine_dimensions_D(deg, fList, time_multiplier, hash_table, opts.use_oldhash);
end

end
