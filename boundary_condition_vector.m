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
nTerms = numel(pde.terms);
nDims = numel(dims);

bcVec = zeros(deg^nDims*num_elements,1);

for tt = 1:nTerms % Construct a BC object for each term
    
    term_nD = terms{tt};
    
    for d1 = 1:nDims
        
        dim = dims{d1};
        term_1D = term_nD.terms_1D{d1};
        %         type = term_1D_.type;
        
        xMin = dim.domainMin;
        xMax = dim.domainMax;
        FMWT = dim.FMWT;
        
        lev = dim.lev;
        %         N_1D = 2^lev;
        %         dof_1D = deg * N_1D;
        
        %%
        % Here we account for the different types of operators which do
        % (grad,diff) or do not (mass) require BCs
        
        %         n_parts = 0;
        %         if strcmp(type,'grad')
        %             n_parts = 1;
        %             n_types = {type};
        %             n_BCLs = {dim.BCL};
        %             n_BCRs = {dim.BCR};
        %             n_BCL_fLists = {dim.BCL_fList};
        %             n_BCR_fLists = {dim.BCR_fList};
        %         elseif strcmp(type,'diff')  % to allow for diff being composed of two grad pieces
        %             n_parts = 2;
        %             n_types = {'grad','grad'};
        %             n_BCLs = {term_1D_.BCL1,term_1D_.BCL2};
        %             n_BCRs = {term_1D_.BCR1,term_1D_.BCR2};
        %             n_BCL_fLists = {term_1D_.BCL1_fList,term_1D_.BCL2_fList};
        %             n_BCR_fLists = {term_1D_.BCR1_fList,term_1D_.BCR2_fList};
        %         end
        
        for p=1:numel(term_1D.pterms)
            
            this_type = term_1D.pterms{p}.type;
            
            if strcmp(this_type,'grad') % BCs are only present for grad terms
                
                this_BCL = term_1D.pterms{p}.BCL;
                this_BCR = term_1D.pterms{p}.BCR;
                BCL_fList = term_1D.pterms{p}.BCL_fList;
                BCR_fList = term_1D.pterms{p}.BCR_fList;
                
                %%
                % Initialize to zero
                
                for d2=1:nDims
                    this_dof_1D = deg * 2^dims{d2}.lev;
                    bcL{d1}{d2} = zeros(this_dof_1D,1);
                    bcR{d1}{d2} = zeros(this_dof_1D,1);
                end
                
                timeFacL = 1;
                timeFacR = 1;
                
                if strcmp(this_type,'grad')
                    
                    if strcmp(this_BCL,'D') % Left side
                        
                        %%
                        % Get time multiplier
                        
                        timeFacL = BCL_fList{nDims+1}(time);
                        
                        %%
                        % Get boundary functions for all dims
                        
                        for d2=1:nDims
                            bcL{d1}{d2} = forward_wavelet_transform(pde.deg,pde.dimensions{d2}.lev,...
                                pde.dimensions{d2}.domainMin,pde.dimensions{d2}.domainMax,...
                                BCR_fList{d2},pde.params,time);
                        end
                        
                        %%
                        % Overwrite the trace (boundary) value just for this dim
                        % Func*v|_xMin and Func*v|_xMax
                        
                        bcL_tmp = compute_boundary_condition(pde,time,lev,deg,xMin,xMax,BCL_fList{d1},'L');
                        bcL_tmp = FMWT * bcL_tmp;
                        
                        %%
                        % LDG requires and additional step
                        
                        %                         if strcmp(type,'diff')
                        %                             bcL_tmp = term_1D.mat1*bcL_tmp;  % TODO : is this mat1 always? or mat2 sometimes?
                        %                         end
                        
                        if p > 1
                            bcL_tmp = term_1D.pterms{p-1}.mat * bcL_tmp;
                        end
                        
                        bcL{d1}{d1} = bcL_tmp;
                    end
                    
                    if strcmp(this_BCR,'D') % Right side
                        
                        %%
                        % Get time multiplier
                        
                        timeFacR = BCR_fList{nDims+1}(time);
                        
                        %%
                        % Get boundary functions for all dims
                        
                        for d2=1:nDims
                            bcR{d1}{d2} = forward_wavelet_transform(pde.deg,pde.dimensions{d2}.lev,...
                                pde.dimensions{d2}.domainMin,pde.dimensions{d2}.domainMax,...
                                BCR_fList{d2},pde.params,time);
                        end
                        
                        %%
                        % Overwrite the trace (boundary) value just for this dim
                        % Func*v|_xMin and Func*v|_xMax
                        
                        bcR_tmp = compute_boundary_condition(pde,time,lev,deg,xMin,xMax,BCR_fList{d1},'R');
                        bcR_tmp = FMWT * bcR_tmp;
                        
                        %%
                        % LDG requires and additional step
                        %
                        %                         if strcmp(type,'diff')
                        %                             bcR_tmp = term_1D.mat1*bcR_tmp; % TODO : is this mat1 always? or mat2 sometimes?
                        %                         end
                        
                        if p > 1
                            bcR_tmp = term_1D.pterms{p-1}.mat * bcR_tmp;
                        end
                        
                        bcR{d1}{d1} = bcR_tmp;
                    end
                    
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