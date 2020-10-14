%% Re-chain Time-Independent 1D coefficient matrices
% If max level coefficient construction is used, partial term coefficients
% are constructed for some maximum adaptivity level. then, when adaptivity
% dictates (enables/disables elements that increase/decrease level for some 
%dimension(s)), we re-chain the stored partial term coefficients for new level(s).

function pde = get_coeff_mats_rechain(pde, deg, new_levels)

num_terms = numel(pde.terms);
num_dims = numel(pde.dimensions);

assert(size(new_levels, 1) == num_dims);


for dim = 1:num_dims
    % if the new level is different than old level, we will need to rechain
    if pde.dimensions{dim}.lev ~= new_levels(dim)
        for term = 1:num_terms
            if ~pde.terms{term}.terms_1D{dim}.time_dependent 
                %TD terms will be regenerated anyway
            [pde.terms{term}.terms_1D{dim}.mat, ...
             pde.terms{term}.terms_1D{dim}.mat_unrotated] = ...
               rechain(pde.terms{term}.terms_1D{dim}, new_levels(dim), deg);
            
            end
        end
        
        if ~isempty(pde.termsLHS)    
            num_terms_left = numel(pde.termsLHS);
            for term = 1:num_terms_left
            if ~pde.termsLHS{term}.terms_1D{dim}.time_dependent.time_dependent
                [pde.termsLHS{term}.terms_1D{dim}.mat, ...
                 pde.termsLHS{term}.terms_1D{dim}.mat_unrotated] = ...
                    rechain(pde.termsLHS{term}.terms_1D{dim}, new_levels(dim), deg);
            end
            end
        end
    end
end


end

function [mat, mat_unrotated] = rechain(term_1D, new_level, degree)

n = degree * 2^new_level;
mat = eye(n);
mat_unrotated = eye(n);

for i=1:numel(term_1D.pterms)
    pterm_mat = term_1D.pterms{i}.mat;
    assert(size(pterm_mat, 1) >= n);
    assert(size(pterm_mat, 2) >= n);
    mat = mat * pterm_mat(1:n, 1:n);
    
    pterm_mat_unrot = term_1D.pterms{i}.mat_unrotated;
    assert(size(pterm_mat_unrot, 1) >= n);
    assert(size(pterm_mat_unrot, 2) >= n);
    mat_unrotated = mat_unrotated * pterm_mat_unrot(1:n, 1:n);
    
end

end
