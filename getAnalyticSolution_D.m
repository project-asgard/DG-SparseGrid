function value = getAnalyticSolution_D(x,t,pde)

%%
% Here x is a D dimensional cell array, i.e,
% for a 2D (x,v) problem, you might do
% ans = getAnalyticSolution({1.0,5.0},10.0) for x=1,v=5,t=10.
% Make sense? Good.

dimensions = pde.dimensions;
a = pde.analytic_solutions_1D;

nDims = numel(dimensions);
assert(numel(x)==nDims,'Input coordinate has incorrect dimensionality');
assert(numel(a)==nDims+1,'Analytic solution specified in PDE spec has incorrect dimensionality');

value = 1;

%%
% Apply 1D coordinate variations

for d=1:nDims
    foo = a{d}; % Function handle
    value = value .* foo(x{d});
end

%%
% Apply time variation
foo = a{nDims+1};
value = value .* foo(t);

end