function bcVec = getBoundaryCondition1(pde,HashInv,time)

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

deg = dims{1}.deg; % TODO

nHash = numel(HashInv);
nTerms = numel(pde.terms);
nDims = numel(dims);

bcVec = zeros(deg^nDims*nHash,1);

for tt = 1:nTerms % Construct a BC object for each term
    
    term = terms{tt};
    
    for d1 = 1:nDims
        
        dim1 = dims{d1};
        
        xMin = dim1.domainMin;
        xMax = dim1.domainMax;
        FMWT = dim1.FMWT;
        BCL = dim1.BCL;
        BCR = dim1.BCR;
        BCL_fList = dim1.BCL_fList;
        BCR_fList = dim1.BCR_fList;
        
        ftL = 1;
        ftR = 1;
        if strcmp(BCL,'D') % dirichlet
            ftL = dim1.BCL_fList{nDims+1}(time);
        end
        if strcmp(BCR,'D') % dirichlet
            ftR = dim1.BCR_fList{nDims+1}(time);
        end
        
        lev = dim1.lev;
        N = 2^lev;
        dof_1D = deg * N;
        
        %%
        % Initialize to zero for when the term is of type 2 for example
        
        for d2=1:nDims
            bcL{d1}{d2} = zeros(dof_1D,1);
            bcR{d1}{d2} = zeros(dof_1D,1);
        end
        
        %%
        % Get the boundary integral for all other dimensions
        
        if strcmp(term{d1}.type,'grad') || strcmp(term{d1}.type,'diff') % grad or diffusion operators
            
            if strcmp(BCL,'D')
                bcL{d1} = ComputeRHS(pde,time,nDims,dim1,'L'); % returns a nDim length list
            end
            if strcmp(BCR,'D')
                bcR{d1} = ComputeRHS(pde,time,nDims,dim1,'R'); % returns a nDim length list
            end
            
        end
        
        
        %%
        % Overwrite the boundary value just for this dimension
        
        % Func*v|_xMin and Func*v|_xMax
        
        if strcmp(term{d1}.type,'grad') || strcmp(term{d1}.type,'diff') % grad or diffusion operators
            
            if strcmp(BCL,'D') % Dirichlet
                
                bcL_tmp = ComputeBC(pde,time,lev,deg,xMin,xMax,BCL_fList{d1},'L');
                bcL_tmp = FMWT * bcL_tmp;
                
                if strcmp(term{d1}.type,'diff') % LDG requires additional step
                    bcL_tmp = term{d1}.matD*bcL_tmp;
                end
                
                bcL{d1}{d1} = bcL_tmp;
            end
            
            if strcmp(BCR,'D') % Dirichlet
                
                bcR_tmp = ComputeBC(pde,time,lev,deg,xMin,xMax,BCR_fList{d1},'R');
                bcR_tmp = FMWT * bcR_tmp;
                
                if strcmp(term{d1}.type,'diff') % LDG requires additional step
                    bcR_tmp = term{d1}.matD*bcR_tmp;
                end
                
                bcR{d1}{d1} = bcR_tmp;
            end
            
        end
        
        fListL = bcL{d1};
        fListR = bcR{d1};
        bcVec = bcVec + combine_dimensions_D(fListL,ftL,HashInv,pde);
        bcVec = bcVec + combine_dimensions_D(fListR,ftR,HashInv,pde);
        
    end % loop over dim1
    
end % loop over terms


end