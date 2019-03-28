function [bc3] = getBoundaryConditionVectors(pde,HashInv,bcL,bcR,time)

nDims = numel(pde.dimensions);
nHash = numel(HashInv);
deg = pde.dimensions{1}.deg;

for d = 1:nDims
     
    dim = pde.dimensions{d};
    
    lev = dim.lev;
    
    bcL1_tmp = ComputeRHS(nDims,lev,deg,dim,dim.BCL_fList);
    bcR1_tmp = ComputeRHS(nDims,lev,deg,dim,dim.BCR_fList);
    
    for d2 = 1:nDims % only compute for the dimension other than d
        if abs(d2-d)>0
            bcL1{d} = (dim.FMWT*bcL1_tmp{d2});
            bcR1{d} = (dim.FMWT*bcR1_tmp{d2});
        end
    end
    
end

%%
% Apply time dependence

for d=1:nDims
    
    dim = pde.dimensions{d};
    
    bcL_t{d} = bcL{d} * dim.BCL_fList{nDims+1}(time);
    bcR_t{d} = bcR{d} * dim.BCR_fList{nDims+1}(time);

    bcL1{d}
    
end

%%
% Construct the RHS piece of the Dirichlet BCs

ft = 1;
bc3 = zeros(deg^nDims*nHash,1);
for d=1:nDims
    
    %%
    % For each dimensional boundary
    clear fListL fListR;
    for d2=1:nDims
        
        if d == d2
            fListL{d2} = bcL{d};
            fListR{d2} = bcR{d};
        else
            fListL{d2} = bcL1{d2};
            fListR{d2} = bcR1{d2};
        end
        
    end
    
    bc3 = bc3 + combine_dimensions_D(fListL,ft,HashInv,pde);
    bc3 = bc3 + combine_dimensions_D(fListR,ft,HashInv,pde);
    
end

%%
% Construct the RHS piece of the Neumann BCs
%
% TODO
%
end