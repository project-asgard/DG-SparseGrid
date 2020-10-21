function term_out = term_fill(term_in)

term_out = term_in;

num_dims = numel(term_in);

for d=1:num_dims
    
    %%
    % Check if this dim for this term is empty. If so, populate with mass
    % matrix.
    
    if isempty(term_in{d})
        
        term_out{d} = term_mass;
        term_out{d} = term_1D({mass()});
        
    end
    
end

end
