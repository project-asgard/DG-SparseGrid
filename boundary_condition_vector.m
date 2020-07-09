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

deg = pde.deg;

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
        %         type = term_1D_.type;
        
        xMin = dim.domainMin;
        xMax = dim.domainMax;
        FMWT = dim.FMWT;
        
        lev = dim.lev;
        
        for p=1:numel(term_1D.pterms)
            
            this_type = term_1D.pterms{p}.type;
            this_g = term_1D.pterms{p}.g;
            
            if strcmp(this_type,'grad') % BCs are only present for grad terms
                
                this_BCL = term_1D.pterms{p}.BCL;
                this_BCR = term_1D.pterms{p}.BCR;
                BCL_fList = term_1D.pterms{p}.BCL_fList;
                BCR_fList = term_1D.pterms{p}.BCR_fList;
                
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
                    
                    timeFacL = BCL_fList{num_dims+1}(time);
                    
                    %%
                    % Get boundary functions for all dims
                    
                    for d2=1:num_dims
                        bcL{d1}{d2} = forward_wavelet_transform(pde.deg,pde.dimensions{d2}.lev,...
                            pde.dimensions{d2}.domainMin,pde.dimensions{d2}.domainMax,...
                            BCL_fList{d2},pde.params,time);
                    end
                    
                    %%
                    % Overwrite the trace (boundary) value just for this dim
                    % Func*v|_xMin and Func*v|_xMax
                    
                    bcL_tmp = compute_boundary_condition(pde,this_g,time,lev,deg,xMin,xMax,BCL_fList{d1},'L');
                    bcL_tmp = FMWT * bcL_tmp;
                    
                    %%
                    % Apply mats from preceeding pterms when chaining (p>1)
                    % FIXME : test this for p>2
                    
                    if p > 1
                        preceeding_mat = eye(size(term_1D.pterms{1}.mat));
                        for nn=1:p-1
                            preceeding_mat = preceeding_mat * term_1D.pterms{nn}.mat;
                        end
                        bcL_tmp = preceeding_mat * bcL_tmp;
                    end
                    
                    bcL{d1}{d1} = bcL_tmp;
                end
                
                if strcmp(this_BCR,'D') % Right side
                    
                    %%
                    % Get time multiplier
                    
                    timeFacR = BCR_fList{num_dims+1}(time);
                    
                    %%
                    % Get boundary functions for all dims
                    
                    for d2=1:num_dims
                        bcR{d1}{d2} = forward_wavelet_transform(pde.deg,pde.dimensions{d2}.lev,...
                            pde.dimensions{d2}.domainMin,pde.dimensions{d2}.domainMax,...
                            BCR_fList{d2},pde.params,time);
                    end
                    
                    %%
                    % Overwrite the trace (boundary) value just for this dim
                    % Func*v|_xMin and Func*v|_xMax
                    
                    bcR_tmp = compute_boundary_condition(pde,this_g,time,lev,deg,xMin,xMax,BCR_fList{d1},'R');
                    bcR_tmp = FMWT * bcR_tmp;
                    
                    %%
                    % Apply mats from preceeding terms when chaining (p>1)
                    
                    if p > 1
                        preceeding_mat = eye(size(term_1D.pterms{1}.mat));
                        for nn=1:p-1
                            preceeding_mat = preceeding_mat * term_1D.pterms{nn}.mat;
                        end
                        bcR_tmp = preceeding_mat * bcR_tmp;
                    end
                    
                    bcR{d1}{d1} = bcR_tmp;
                end
                                
                fListL = bcL{d1};
                fListR = bcR{d1};
                
                bcVec = bcVec + combine_dimensions_D(pde.deg,fListL,timeFacL,hash_table,opts.use_oldhash);
                bcVec = bcVec + combine_dimensions_D(pde.deg,fListR,timeFacR,hash_table,opts.use_oldhash);
                
            end
            
        end
        
    end % loop over dim1
    
end % loop over terms


end