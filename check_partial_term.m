function pterm_out = check_partial_term(num_dimensions,pterm)

% Check to make sure each partial term has the all the right fields.

default_term_mass
default_term_grad
default_term_diff

for t=1:numel(pterm)
    
    this_pterm = pterm{t};
    
    if strcmp(this_pterm.type,'mass')
        this_pterm_out = term_mass;
    end
    if strcmp(this_pterm.type,'grad')
        this_pterm_out = term_grad;
    end
    if strcmp(this_pterm.type,'diff')
        this_pterm_out = term_diff;
    end
    
    this_pterm_out.coeff_mat = [];
    this_pterm_out.coeff_mat0 = [];
    this_pterm_out.mat1 = [];
    this_pterm_out.mat2 = [];
    
    % Check to make sure all fields exist.
    % If not, use default.
    
    fn = fieldnames(this_pterm_out);
    for k=1:numel(fn)
        if isfield(this_pterm,fn{k})
            this_pterm_out.(fn{k}) = this_pterm.(fn{k});
        end
    end
    
    %%
    % Check if there are erroneous field names
    
    fn = fieldnames(this_pterm);
    for k=1:numel(fn)
        if ~isfield(this_pterm_out,fn{k})
            error(strcat('Unrecognized term in term: ', fn{k} ));
        end
    end
    
    pterm_out{t} = this_pterm_out;
    
end

end