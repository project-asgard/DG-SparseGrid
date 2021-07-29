function bcVec = boundary_condition_vector(pde, opts, hash_table, time)

%%
% If inhomogeneous Dirichlet, add a combined left and right term which
% is addative, i.e., if only left is D then the right component is added in
% as zero.

%%
% pde shortcuts

dims = pde.dimensions;
terms = pde.terms;

%%
% dim shortcuts

deg = opts.deg;

if opts.use_oldhash
    num_elements = numel(hash_table);
else
    num_elements = numel(hash_table.elements_idx);
end
num_terms = numel(pde.terms);
num_dims = numel(dims);

bcVec = zeros(deg^num_dims*num_elements,1);

for tt = 1:num_terms % Construct a BC object for each term
    
    term_nD = terms{tt};
    
    for d1 = 1:num_dims
        
        dim = dims{d1};
        term_1D = term_nD.terms_1D{d1};
        
        xMin = dim.min;
        xMax = dim.max;
        
        lev = dim.lev;
        dof = deg * 2^lev;
        
        for p=1:numel(term_1D.pterms)
            
            pterm = term_1D.pterms{p};
            this_type = pterm.type;
            this_g = pterm.g;
            
            if strcmp(this_type,'grad') || strcmp(this_type,'div') % BCs are only present for grad/div terms
                
                this_BCL = pterm.BCL;
                this_BCR = pterm.BCR;
                BCL_fList = pterm.BCL_fList;
                BCR_fList = pterm.BCR_fList;
                
                %%
                % Initialize to zero
                
                for d2=1:num_dims
                    this_dof_1D = deg * 2^dims{d2}.lev;
                    bcL{d1}{d2} = zeros(this_dof_1D,1);
                    bcR{d1}{d2} = zeros(this_dof_1D,1);
                end
                
                timeFacL = 1;
                timeFacR = 1;
                                
                if strcmp(this_BCL,'D') % Left side
                    
                    %%
                    % Get time multiplier
                    
                    timeFacL = BCL_fList{num_dims+1}(time,pde.params);
                    
                    %%
                    % Get boundary functions for all dims
                    
                    for d2=1:num_dims
                        term_1D_oth_dim = term_nD.terms_1D{d2};
                        mass_pterm      = term_1D_oth_dim.pterms{p};
                        bcL{d1}{d2} = forward_wavelet_transform(opts.deg,pde.dimensions{d2}.lev,...
                            pde.dimensions{d2}.min,pde.dimensions{d2}.max,...
                            BCL_fList{d2},mass_pterm.dV,pde.params,pde.transform_blocks,time);
                        %Apply inverse mat
                        bcL{d1}{d2}     = mass_pterm.LHS_mass_mat \ bcL{d1}{d2};
                        %Apply previous pterms
                        for q=1:p-1
                             bcL{d1}{d2} = term_1D_oth_dim.pterms{q}.mat * bcL{d1}{d2};
                        end
                    end
                    
                    
                    
                    
                    %%
                    % Overwrite the trace (boundary) value just for this dim
                    % Func*v|_xMin and Func*v|_xMax
                    
                    bcL_tmp = compute_boundary_condition(pde,this_g,pterm.dV,time,lev,deg,xMin,xMax,BCL_fList{d1},'L');
                    trans_side = 'LN';
                    bcL_tmp = apply_FMWT_blocks(lev, pde.transform_blocks, bcL_tmp, trans_side);
                    
                    %%
                    % Apply LHS_mass_mat for this pterm
                   
                    M = pterm.LHS_mass_mat;
                    bcL_tmp = M \ bcL_tmp;
                    
                    %%
                    % Apply mats from preceeding pterms when chaining (p>1)
                    
                    preceeding_mat = eye(dof);
                    for nn=1:p-1
                        preceeding_mat = preceeding_mat * term_1D.pterms{nn}.mat;
                    end
                    bcL_tmp = preceeding_mat * bcL_tmp;
                    
                    bcL{d1}{d1} = bcL_tmp;
                end
                
                
                
                if strcmp(this_BCR,'D') % Right side
                    
                    %%
                    % Get time multiplier
                    
                    timeFacR = BCR_fList{num_dims+1}(time,pde.params);
                    
                    %%
                    % Get boundary functions for all dims
                    
                    for d2=1:num_dims
                        term_1D_oth_dim = term_nD.terms_1D{d2};
                        mass_pterm      = term_1D_oth_dim.pterms{p};
                        bcR{d1}{d2} = forward_wavelet_transform(opts.deg,pde.dimensions{d2}.lev,...
                            pde.dimensions{d2}.min,pde.dimensions{d2}.max,...
                            BCR_fList{d2},mass_pterm.dV,pde.params,pde.transform_blocks,time);
                        %Apply inverse mat
                        bcR{d1}{d2}     = mass_pterm.LHS_mass_mat \ bcR{d1}{d2};
                        %Apply previous pterms
                        for q=1:p-1
                             bcR{d1}{d2} = term_1D_oth_dim.pterms{q}.mat * bcR{d1}{d2};
                        end
                    end
                    
                    %%
                    % Overwrite the trace (boundary) value just for this dim
                    % Func*v|_xMin and Func*v|_xMax
                    
                    bcR_tmp = compute_boundary_condition(pde,this_g,pterm.dV,time,lev,deg,xMin,xMax,BCR_fList{d1},'R');
                    trans_side = 'LN';
                    bcR_tmp = apply_FMWT_blocks(lev, pde.transform_blocks, bcR_tmp, trans_side);
                    
                    %%
                    % Apply LHS_mass_mat for this pterm
                   
                    M = pterm.LHS_mass_mat;
                    bcR_tmp = M \ bcR_tmp;                
                    
                    %%
                    % Apply mats from preceeding terms when chaining (p>1)
                                        
                    preceeding_mat = eye(dof);
                    for nn=1:p-1
                        preceeding_mat = preceeding_mat * term_1D.pterms{nn}.mat;
                    end
                    bcR_tmp = preceeding_mat * bcR_tmp;
                    
                    bcR{d1}{d1} = bcR_tmp;
                end
                
                fListL = bcL{d1};
                fListR = bcR{d1};
                
                bcVec = bcVec + combine_dimensions_D(opts.deg,fListL,timeFacL,hash_table,opts.use_oldhash);
                bcVec = bcVec + combine_dimensions_D(opts.deg,fListR,timeFacR,hash_table,opts.use_oldhash);
                
            end
            
        end
        
    end % loop over dim1
    
end % loop over terms


end
