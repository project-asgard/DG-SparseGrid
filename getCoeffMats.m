function pde = getCoeffMats (pde, t, TD)

%%
% t : time passed to G function
% TD == 1 for time dependent
% TD == 0 for time independent

nTerms = numel(pde.terms);
nDims = numel(pde.dimensions);

if TD
    TD_STR = 'TD';
else
    TD_STR = 'TI';
end

for tt = 1:nTerms
    
    term = pde.terms{tt};
    
    %%
    % Add dim many operator matrices to each term.
    for d = 1:nDims
        
        dim = pde.dimensions{d};
        
        if term{d}.TD == TD
            
            disp([TD_STR ' - term : ' num2str(tt) '  d : ' num2str(d) ]);
            
            [mat,matD] = coeff_matrix2(pde,t,dim,term{d});
            
            pde.terms{tt}{d}.coeff_mat = mat;
            if strcmp(term{d}.type,'diff') % Keep matU and matD from LDG for use in BC application
                pde.terms{tt}{d}.matD = matD;    
            end
            
        end
    end
end

end