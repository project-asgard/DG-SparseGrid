function [pde] = compute_dimension_mass_mat(opts,pde)
%Computing mass matrix using the moment_dV volume jacobian.
%Needed when determining L2 projections with non-cartesian coordinates.  

num_dims = numel(pde.dimensions);
dims = pde.dimensions;


for d = 1:num_dims
    dim = dims{d};
    construction_level = dim.lev;
    if opts.max_lev_coeffs && ~term_1D.time_dependent
        construction_level = opts.max_lev;
    end
    g1 = @(x,p,t,d) 0*x+1;
    pterm = MASS(g1,[],[],dim.moment_dV);
    [M,~] = coeff_matrix(opts.deg,0,dim,pterm,pde.params,pde.transform_blocks,construction_level);
    pde.dimensions{d}.mass_mat = M;
end

end

