function [pde,bcL,bcR] = getCoeffMats (pde, t, TD)

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

dims = pde.dimensions;

for d=1:nDims
    lev = dims{d}.lev;
    deg = dims{d}.deg;
    N = 2^lev;
    dof_1D = deg * N;
    bcL{d} = zeros(dof_1D,1);
    bcR{d} = zeros(dof_1D,1);
end

for term = 1:nTerms
    
    thisTerm = pde.terms{term};
    
    %%
    % Add dim many operator matrices to each term.
    for d = 1:nDims
        
        if thisTerm{d}.TD == TD
            
            disp([TD_STR ' - term : ' num2str(term) '  d : ' num2str(d) ]);
            
            %coeff_mat = coeff_matrix(t,pde.dimensions{d},pde.terms{term}{d});
            
            [coeff_mat,bcL_tmp,bcR_tmp] = coeff_matrix2(t,pde.dimensions{d},pde.terms{term}{d},d);
            bcL{d} = bcL{d} + bcL_tmp; % TODO : are these additive? zero valued for mass matrix?
            bcR{d} = bcR{d} + bcR_tmp;
            
            pde.terms{term}{d}.coeff_mat = coeff_mat;
            
        end
    end
end

end