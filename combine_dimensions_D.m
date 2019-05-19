function [fval] = combine_dimensions_D(fD,ft,HASHInv,pde)

% Combine (via kron product) a set of 1D multiwavelet transforms to form
% the higher D sparse-grid multiwavelet representation.

% ft is the time multiplier.

nDims = numel(fD);
if pde.useHash
    N = numel(HASHInv);
else
    N = numel(pde.elementsIDX);
end

Deg = pde.dimensions{1}.deg; % TODO Need to generalize this deg_D

fval = sparse(Deg^nDims * N,1);

% nze = nonzeros(pde.elements.coords{1}.lev);

for i=1:N
    
    %%
    % Kron product approach
    
    clear kronMatList;
    for d=1:nDims
        
        if pde.useHash
            ll=HASHInv{i};
            ID = ll(nDims*2+d); % TODO : Check if this indexing correct for D != 2?
        else
            IDlev = pde.elements.lev_p1(pde.elementsIDX(i),d)-1;
            IDpos = pde.elements.pos_p1(pde.elementsIDX(i),d)-1;
            IDe = lev_cell_to_singleD_index(IDlev,IDpos);
            ID = IDe;
        end
        
        index_D = [(ID-1)*Deg+1 : ID*Deg];
        f = fD{d};
        ftmp = f(index_D);
        kronMatList{d} = ftmp;
    end
    
    B = ft;
    
    use_krond = 0;
    if use_krond
        A = krond(nDims,kronMatList);
    else
        A = 1;
        for d=1:nDims
            A = kron(A,kronMatList{d});
        end
    end
    
    tmp = A * B;
    
    Index = Deg^nDims*(i-1)+1:Deg^nDims*i;
    fval(Index,1) = fval(Index,1) + tmp(:);
    
end

end


