function total_value = get_analytic_realspace_solution_D(pde,opts,x,t)

%%
% Here x is a D dimensional cell array, i.e,
% for a 2D (x,v) problem, you might do
% ans = getAnalyticSolution({1.0,5.0},10.0) for x=1,v=5,t=10.
% Make sense? Good.

dims = pde.dimensions;

num_dims = numel(dims);
assert(numel(x)==num_dims,...
    'Input coordinate has incorrect dimensionality');


%%
% Apply 1D coordinate variations

if opts.many_solution_capable
    
    num_solutions = numel(pde.solutions);
    
    total_value = 0;
    
    for s = 1:num_solutions
        
        value = 1;
        
        for d=1:num_dims
            foo = pde.solutions{s}{d}; % Function handle
            value = value .* foo(x{d},pde.params,t);
        end
        
        %%
        % Apply time variation
        foo = pde.solutions{s}{num_dims+1};
        value = value .* foo(t);
        
        total_value = total_value + value;
        
    end
    
else
    
    a = pde.analytic_solutions_1D;   
    assert(numel(a)==num_dims+1,...
        'Analytic solution specified in PDE spec has incorrect dimensionality');
    
    value = 1;
    
    for d=1:num_dims
        foo = a{d}; % Function handle
        value = value .* foo(x{d},pde.params,t);
    end
    
    %%
    % Apply time variation
    foo = a{num_dims+1};
    value = value .* foo(t);
    
    total_value = value;
    
end

end