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
    N = numel(HASHInv);
    Ne = numel(pde.elementsIDX);
    assert(N==Ne);
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
fnew = sparse(dof_1D_FG^nDims,1);
f_rSpace = sparse(dof_1D_FG^nDims,1);

%% 
% To be removed
if nDims==2
    A = Meval_D{1};
    B = Meval_D{2};
    f = f_wSpace;
end

for i=1:N
    
    ll=HASHInv{i};
    
    %%
    % Retain 2D version for checking but don't use
    if nDims==2
        I1=HASHInv{i}(5);
        I2=HASHInv{i}(6);
        
        index_I1=[(I1-1)*deg+1:I1*deg];
        index_I2=[(I2-1)*deg+1:I2*deg];
        
        Index = deg^2*(i-1)+1:deg^2*i;
        
        tmp=kron(...
            A(:,index_I1),...
            B(:,index_I2) ...
            )*f(Index(:));
        fnew=fnew+tmp;
    end
    
    %%
    % kron(tmp1,tmp2,...,tmpD)*fwav
    
    clear kronMatList;
    for d=1:nDims
        
        if pde.useHash
            ID = ll(nDims*2+d);
        else
            ID = ll(nDims*2+d);
            %% Etable
            IDlev  = pde.elements{d}.lev(pde.elementsIDX(i));
            IDcell = pde.elements{d}.cell(pde.elementsIDX(i));
            IDe = lev_cell_to_singleD_index(IDlev-1,IDcell-1);
            assert(ID==IDe);
            ID = IDe;
        end
        index_D = [(ID-1)*deg+1 : ID*deg];
       
        %%
        % Assert match the 2D case
        if nDims==2
            if d==1
                assert(norm(index_D-index_I1)==0);
            end
            if d==2
                assert(norm(index_D-index_I2)==0);
            end
        end
        
        thisMat = Meval_D{d};
        thisMat1 = thisMat(:,index_D);
        
        %%
        % Assert match the 2D case
        if nDims==2
            if d==1
                assert(norm(full(thisMat1-A(:,index_I1)))==0);
            end
            if d==2
                assert(norm(full(thisMat1-B(:,index_I2)))==0);
            end
        end
        
        kronMatList{d} = thisMat1; % Build kron list
    end
   
    element_ii = deg^nDims*(i-1)+1:deg^nDims*i;
    
    X = f_wSpace(element_ii);

    if use_kronmultd
        Y = kron_multd(nDims,kronMatList,X);
    else
        Y = kron_multd_full(nDims,kronMatList,X);
    end
        
    %%
    % Test the kron_multd
    
    AA = 1;
    for d=1:nDims
        AA = kron(AA,kronMatList{d});  
    end
    YA = AA * X;
    
    tol=1e-12;
    assert(norm(YA-Y)<tol);
    
    tol = 1e-12;
    if nDims==2
        assert(norm(Y-tmp)<tol);
    end
    
    f_rSpace = f_rSpace + Y;
    
end

if nDims==2
    assert(norm(fnew-f_rSpace)<tol);
end

end



