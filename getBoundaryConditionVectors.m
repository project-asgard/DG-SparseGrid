function [bc3] = getBoundaryConditionVectors(pde,HashInv,bcL,bcR,time)

nDims = numel(pde.dimensions);
nHash = numel(HashInv);
deg = pde.dimensions{1}.deg;

% time = 0;
for d = 1:nDims
    
%     [Delta,bcL_tmp,bcR_tmp] = coeff_matrix2(t,pde.dimensions{d},pde.terms{1}{1},d);
%     bcL{d} = bcL_tmp;
%     bcR{d} = bcR_tmp;
    
    dim = pde.dimensions{d};
    
    lev = dim.lev;
    
    bcL1_tmp = ComputeRHS(lev,deg,dim,dim.BCL_fList,time);
    bcR1_tmp = ComputeRHS(lev,deg,dim,dim.BCR_fList,time);
    
    for d2 = 1:nDims % only compute for the dimension other than d
        if abs(d2-d)>0
            bcL1{d} = (dim.FMWT*bcL1_tmp{d2});
            bcR1{d} = (dim.FMWT*bcR1_tmp{d2});
        end
    end
    
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