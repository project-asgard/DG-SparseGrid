function fval = md_eval_function(opts, deg, dims, params, md_funcs, hash_table, blocks, t, mass_bool)

% Returns a multi-D wavelet space function defined by "md_funcs"

% This is done by a standard L^2 projection: fval solves the system
%
% M*fval = RHS
%
% where RHS is the vector defined by RHS_i = (md_funcs,phi_i)_\W 
%
% If the coordinate system is cartesian, then M=I and fval=RHS.  Else we
% need to invert by the mass matrix, but since M can be written as a
% kronecker product, we can use instead invert the mass matrix 
% on each dimension

% mass_bool = 1 implies that fval is returned.
% mass_bool = 0 implies that M*fval is returned.  No mass matrix inversion
% is done -- this is useful when calculating moment integrals for
% conservation property.  

if nargin < 9 %mass_bool not added
    mass_bool = 1;
end

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
            md_func{d},dims{d}.moment_dV,params, blocks, t);
        if mass_bool
            N = numel(fList{d});
            fList{d} = dims{d}.mass_mat(1:N,1:N) \ fList{d}; 
        end
    end
    
    time_multiplier = md_funcs{s}{num_dims+1}(t);
    fval = fval + combine_dimensions_D(deg, fList, time_multiplier, hash_table, opts.use_oldhash);
end

end
