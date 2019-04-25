function pde = checkTerms(pde)

terms = pde.terms;
dims = pde.dimensions;

ndims = numel(dims);
nterms = numel(terms);

for t=1:nterms
    for d=1:ndims
        terms{t}{d} = checkPartialTerm(ndims,terms{t}{d});
    end
end

pde.terms = terms;

end