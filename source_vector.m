function fval = source_vector(pde,opts,hash_table,t)

% Returns the wavelet transformed source


% num_dimensions  = numel(pde.dimensions);
% num_sources     = numel(pde.sources);
% 
% %%
% % Loop over the number of sources, each of which has nDims + time elements.
% 
% fval = 0;
% for s=1:num_sources
%     for d=1:num_dimensions
%         fList{d} = forward_wavelet_transform(opts.deg,pde.dimensions{d}.lev,...
%             pde.dimensions{d}.min,pde.dimensions{d}.max,...
%             pde.sources{s}{d},pde.params,time);
%     end
%     
%     time_multiplier = pde.sources{s}{num_dimensions+1}(time);
%     fval = fval + combine_dimensions_D(opts.deg, fList, time_multiplier, hash_table, opts.use_oldhash);
% end

fval = md_eval_function(opts, opts.deg, pde.dimensions, ...
    pde.params, pde.sources, hash_table, pde.transform_blocks, t);

end
