function pde = get_coeff_mats (pde,opts,t,TD,use_oldcoeffmat)

%%
% t : time passed to G function
% TD == 1 for time dependent
% TD == 0 for time independent

num_terms = numel(pde.terms);
num_dimensions = numel(pde.dimensions);

oldcoeff = 0;
if exist('use_oldcoeffmat','var') && ~isempty(use_oldcoeffmat)
    oldcoeff = use_oldcoeffmat;
end

debug = 0;

if TD
    TD_STR = 'TD';
else
    TD_STR = 'TI';
end

%%
% Normal RHS terms

for tt = 1:num_terms
    
    term_nD = pde.terms{tt};
    
    %%
    % Add dim many operator matrices to each term.
    for d = 1:num_dimensions
        
        dim = pde.dimensions{d};
        term_1D = term_nD.terms_1D{d};

        if term_1D.time_dependent == TD
            
            if debug; disp([TD_STR ' - term : ' num2str(tt) '  d : ' num2str(d) ]); end           
            
            if oldcoeff
                mat = coeff_matrix_old(opts.deg,t,dim,term_1D);
                pde.terms{tt}{d}.coeff_mat = mat;
            else
                construction_level = dim.lev;
                if opts.max_lev_coeffs && ~term_1D.time_dependent
                    construction_level = opts.max_lev;
                end
                    
                [term_1D_out] = coeff_matrix(num_dimensions,opts.deg,t,dim,term_1D,pde.params, ...
                                             pde.transform_blocks, construction_level);
                pde.terms{tt}.terms_1D{d} = term_1D_out;

            end
            
            
        end
    end
end

%%
% LHS mass matrix

if ~isempty(pde.termsLHS)
    
    nTermsLHS = numel(pde.termsLHS);
    
    for tt=1:nTermsLHS
        
        term_nD = pde.termsLHS{tt};
        
        for d = 1:num_dimensions
            
            dim = pde.dimensions{d};
            term_1D = term_nD.terms_1D{d};
            
            if term_1D.time_dependent == TD
                
                if debug; disp([TD_STR ' - LHS term : ' num2str(1) '  d : ' num2str(d) ]); end
                
                for p=1:numel(term_1D.pterms)
                    assert(strcmp(term_1D.pterms{p}.type,'mass'));
                end
                
                if oldcoeff
                    error('Non-identity LHS mass matrix not supported by "use_oldcoeffmat=1"');
                else          
                    construction_level = dim.lev;
                    if opts.max_lev_coeffs && ~term_1D.time_dependent
                        construction_level = opts.max_lev;
                    end
                    [term_1D_out] = coeff_matrix(num_dimensions,opts.deg,t,dim,term_1D,pde.params,...
                        pde.transform_blocks, construction_level);
                    
                    pde.termsLHS{tt}.terms_1D{d} = term_1D_out;
                end
                
            end
            
        end
    end
    
end

end
