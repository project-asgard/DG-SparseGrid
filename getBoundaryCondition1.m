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

deg = dims{1}.deg; % TODO : assumes deg independent of dim

if pde.useHash
    N = numel(HashInv);
else
%     N = numel(HashInv);
    Ne = numel(pde.elementsIDX);
%     assert(N==Ne);
    N=Ne;
end
nTerms = numel(pde.terms);
nDims = numel(dims);

bcVec = zeros(deg^nDims*N,1);

for tt = 1:nTerms % Construct a BC object for each term
    
    term = terms{tt};
    
    for d1 = 1:nDims
        
        dim = dims{d1};
        pTerm = term{d1};
        type = pTerm.type;
        
        xMin = dim.domainMin;
        xMax = dim.domainMax;
        FMWT = dim.FMWT;
        
%         lev = dim.lev;
        lev = max(pde.elements.lev_p1(:,d1)-1);
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
                bcL{d1}{d2} = zeros(dof_1D,1);
                bcR{d1}{d2} = zeros(dof_1D,1);
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
                    
%                     bcL{d1} = ComputeRHS(pde,time,dim,BCL_fList,FMWT); % returns a nDim length list
                    for d2=1:nDims
                        bcL{d1}{d2} = forwardMWT(pde,d2,BCR_fList{d2},time);
                    end
                    
                    %%
                    % Overwrite the trace (boundary) value just for this dim
                    % Func*v|_xMin and Func*v|_xMax
                    
                    bcL_tmp = ComputeBC(pde,time,lev,deg,xMin,xMax,BCL_fList{d1},'L');
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
                    
%                     bcR{d1} = ComputeRHS(pde,time,dim,BCR_fList,FMWT); % returns a nDim length list 
                    for d2=1:nDims
                        bcR{d1}{d2} = forwardMWT(pde,d2,BCR_fList{d2},time);
                    end                    
                    
                    %%
                    % Overwrite the trace (boundary) value just for this dim
                    % Func*v|_xMin and Func*v|_xMax
                    
                    bcR_tmp = ComputeBC(pde,time,lev,deg,xMin,xMax,BCR_fList{d1},'R');
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
             
            bcVec = bcVec + combine_dimensions_D(fListL,timeFacL,HashInv,pde);
            bcVec = bcVec + combine_dimensions_D(fListR,timeFacR,HashInv,pde);
            
        end
        
    end % loop over dim1
    
end % loop over terms


end