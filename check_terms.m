function pde = checkTerms(pde)

terms = pde.terms;
termsLHS = pde.termsLHS;
dims = pde.dimensions;

num_dimensions = numel(dims);
nterms = numel(terms);

for t=1:nterms
    for d=1:num_dimensions
        terms{t}{d} = check_partial_term(num_dimensions,terms{t}{d});
    end
end

pde.terms = terms;

ntermsLHS = numel(termsLHS);

for t=1:ntermsLHS
    for d=1:num_dimensions
        termsLHS{t}{d} = check_partial_term(num_dimensions,termsLHS{t}{d});
    end
end

pde.termsLHS = termsLHS;

end