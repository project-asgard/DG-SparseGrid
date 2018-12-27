function [pde] = getCoeffMats (pde, t, TD)

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

    for term = 1:nTerms
        
        thisTerm = pde.terms{term};
        
        %%
        % Add dim many operator matrices to each term.
        for d = 1:nDims
            
            if thisTerm{d}.TD == TD
                
                disp([TD_STR ' - term : ' num2str(term) '  d : ' num2str(d) ]);
                
                coeff_mat = coeff_matrix(t,pde.dimensions{d},pde.terms{term}{d});             
                pde.terms{term}{d}.coeff_mat = coeff_mat;
                
            end
        end
    end
        
end