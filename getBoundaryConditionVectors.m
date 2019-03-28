function [bc3] = getBoundaryConditionVectors(pde,bc,HashInv,time)

nDims = numel(pde.dimensions);
nHash = numel(HashInv);
deg = pde.dimensions{1}.deg;

bcL = bc.bcL;
bcR = bc.bcR;
bcL1 = bc.bcL1;
bcR1 = bc.bcR1;

%%
% Construct the RHS piece of the Dirichlet BCs

bc3 = zeros(deg^nDims*nHash,1);
for d=1:nDims
    
    dim = pde.dimensions{d};
    
    ftL = dim.BCL_fList{nDims+1}(time);
    ftR = dim.BCR_fList{nDims+1}(time);
    
    %%
    % For each dimensional boundary
    clear fListL fListR;
    for d2=1:nDims
        
        if d == d2
            fListL{d2} = bcL{d};
            fListR{d2} = bcR{d};
        else
            fListL{d2} = bcL1{d}{d2};
            fListR{d2} = bcR1{d}{d2};
        end
        
    end
    
    bc3 = bc3 + combine_dimensions_D(fListL,ftL,HashInv,pde);
    bc3 = bc3 + combine_dimensions_D(fListR,ftR,HashInv,pde);
    
end

%%
% Construct the RHS piece of the Neumann BCs
%
% TODO
%
end