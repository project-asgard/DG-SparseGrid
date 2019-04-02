function bcVec = getBoundaryCondition1(pde,HashInv,time)

%%
% If inhomogeneous Dirichlet, add a combined left and right term which
% is addative, i.e., if only left is D then the right component is added in
% as zero.

nHash = numel(HashInv);
nTerms = numel(pde.terms);
nDims = numel(pde.dimensions);
dims = pde.dimensions;
terms = pde.terms;

deg = dims{1}.deg; % TODO

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
        if BCL == 1 % dirichlet
            ftL = dim1.BCL_fList{nDims+1}(time);
        end
        if BCR == 1 % dirichlet
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
        
        if term{d1}.type == 1 || term{d1}.type == 3 % grad or del^2 operators
            
            if BCL == 1
                bcL{d1} = ComputeRHS(nDims,dim1,'L'); % returns a nDim length list
            end
            if BCR == 1
                bcR{d1} = ComputeRHS(nDims,dim1,'R'); % returns a nDim length list
            end
            
        end          

        
        %%
        % Overwrite the boundary value just for this dimension
       
        % Func*v|_xMin and Func*v|_xMax      
        
        if term{d1}.type == 1 || term{d1}.type == 3 % grad or del^2 operators
            
            if BCL == 1 % Dirichlet
                bcL_tmp = ComputeBC(pde,time,lev,deg,xMin,xMax,BCL_fList{d1},'L');
                if term{d1}.type == 3 % LDG requires additional step
                    bcL_tmp = term{d1}.matD*bcL_tmp;
                end
                bcL{d1}{d1} = FMWT * bcL_tmp;
            end
            
            if BCR == 1 % Dirichlet
                bcR_tmp = ComputeBC(pde,time,lev,deg,xMin,xMax,BCR_fList{d1},'R');
                if term{d1}.type == 3 % LDG requires additional step
                    bcR_tmp = term{d1}.matD*bcL_tmp;
                end
                bcR{d1}{d1} = FMWT * bcR_tmp;
            end
            
        end
        
        fListL = bcL{d1};
        fListR = bcR{d1};
        bcVec = bcVec + combine_dimensions_D(fListL,ftL,HashInv,pde);
        bcVec = bcVec + combine_dimensions_D(fListR,ftR,HashInv,pde);
        
    end % loop over dim1
    
end % loop over terms


end