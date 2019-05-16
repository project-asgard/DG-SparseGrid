function [f_rSpace]=Multi_2D_D(pde,Meval_D,f_wSpace,HASHInv)

%%
% Apply dimension many inverse wavelet transforms to the solution vector
% f_wSpace -> f_rSpace (wavelet space -> real space)
% via Kron(tmp1,tmp2,...,tmpD)*f
% Note: the output grid of this transform is specified in matrix_plot_D.m

use_kronmultd = 1;

dimensions = pde.dimensions;

nDims = numel(dimensions);
if pde.useHash
    N = numel(HASHInv);
else
    Ne = numel(pde.elementsIDX);
    N = Ne;
end

deg = dimensions{1}.deg; % TODO : generalize to deg_D
lev = dimensions{1}.lev; % TODO : generalize to lev_D

%%
% Catch for TODO

if nDims>1
    for d=2:nDims
        assert(dimensions{d}.lev==lev);
        assert(dimensions{d}.deg==deg);
    end
end

%%
% TODO : surely this size depends on the parameters in the matrix_plot_D
% routine? i.e., the grid upon which we decide to evaluate the solution?
dof_1D_FG = deg*2^(lev); 

%%
% The real space vector is on the full-grid size? Where is this decided?
% (matrix_plot_D.m I think)
% fnew = sparse(dof_1D_FG^nDims,1);
f_rSpace = sparse(dof_1D_FG^nDims,1);

for i=1:N
    
    if pde.useHash
        ll=HASHInv{i};
    else
    end
    
    clear kronMatList;
    for d=1:nDims
        
        if pde.useHash
            ID = ll(nDims*2+d);
        else
            IDlev = pde.elements.lev(pde.elementsIDX(i),d)-1;
            IDpos = pde.elements.pos(pde.elementsIDX(i),d)-1;
            IDe = lev_cell_to_singleD_index(IDlev,IDpos);
            ID = IDe;
        end
        index_D = [(ID-1)*deg+1 : ID*deg];
        
        thisMat = Meval_D{d};
        thisMat1 = thisMat(:,index_D);
        
        kronMatList{d} = thisMat1; % Build kron list
    end
   
    element_ii = deg^nDims*(i-1)+1:deg^nDims*i;
    
    X = f_wSpace(element_ii);

    if use_kronmultd
        Y = kron_multd(nDims,kronMatList,X);
    else
        Y = kron_multd_full(nDims,kronMatList,X);
    end
        
    f_rSpace = f_rSpace + Y;
    
end

end



