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
    
    term = terms{tt};
    
    for d1 = 1:nDims
        
        dim = dims{d1};
        pTerm = term{d1};
        type = pTerm.type;
        
        xMin = dim.domainMin;
        xMax = dim.domainMax;
        FMWT = dim.FMWT;
        
        lev = dim.lev;
        N_1D = 2^lev;
        dof_1D = deg * N_1D;
        
        %%
        % Here we account for the different types of operators which do
        % (grad,diff) or do not (mass) require BCs
        
        n_parts = 0;
        if strcmp(type,'grad')
            n_parts = 1;
            n_types = {type};
            n_BCLs = {dim.BCL};
            n_BCRs = {dim.BCR};
            n_BCL_fLists = {dim.BCL_fList};
            n_BCR_fLists = {dim.BCR_fList};
        elseif strcmp(type,'diff')  % to allow for diff being composed of two grad pieces
            n_parts = 2;
            n_types = {'grad','grad'};
            n_BCLs = {pTerm.BCL1,pTerm.BCL2};
            n_BCRs = {pTerm.BCR1,pTerm.BCR2};
            n_BCL_fLists = {pTerm.BCL1_fList,pTerm.BCL2_fList};
            n_BCR_fLists = {pTerm.BCR1_fList,pTerm.BCR2_fList};
        end
        
        for n=1:n_parts
            
            thisType = n_types{n};
            thisBCL = n_BCLs{n};
            thisBCR = n_BCRs{n};
            BCL_fList = n_BCL_fLists{n};
            BCR_fList = n_BCR_fLists{n};
            
            %%
            % Initialize to zero
            
            for d2=1:nDims           
                this_dof_1D = deg * 2^dims{d2}.lev;
                bcL{d1}{d2} = zeros(this_dof_1D,1);
                bcR{d1}{d2} = zeros(this_dof_1D,1);
            end
            
            timeFacL = 1;
            timeFacR = 1;
            
            if strcmp(thisType,'grad')
                
                if strcmp(thisBCL,'D') % Left side
                    
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
                    
                    if strcmp(type,'diff')
                        bcL_tmp = pTerm.mat1*bcL_tmp;  % TODO : is this mat1 always? or mat2 sometimes?
                    end
                    
                    bcL{d1}{d1} = bcL_tmp;
                end
                
                if strcmp(thisBCR,'D') % Right side
                    
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
                    
                    if strcmp(type,'diff')
                        bcR_tmp = pTerm.mat1*bcR_tmp; % TODO : is this mat1 always? or mat2 sometimes?
                    end
                    
                    bcR{d1}{d1} = bcR_tmp;
                end
                
            end
            
            fListL = bcL{d1};
            fListR = bcR{d1};
             
            bcVec = bcVec + combine_dimensions_D(pde.deg,fListL,timeFacL,hash_table,opts.use_oldhash);
            bcVec = bcVec + combine_dimensions_D(pde.deg,fListR,timeFacR,hash_table,opts.use_oldhash);
            
        end
        
    end % loop over dim1
    
end % loop over terms


end