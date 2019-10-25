function [f_rSpace]=wavelet_to_realspace(pde, opts, Meval_D, f_wSpace, hash_table)

%%
% Apply dimension many inverse wavelet transforms to the solution vector
% f_wSpace -> f_rSpace (wavelet space -> real space)
% via Kron(tmp1,tmp2,...,tmpD)*f
% Note: the output grid of this transform is specified in matrix_plot_D.m

use_kronmultd = 1;

dimensions = pde.dimensions;

nDims = numel(dimensions);
if opts.use_oldhash
    num_elements = numel(hash_table);
else
    num_elements = numel(hash_table.elements_idx);
end

deg = pde.deg; % TODO : generalize to deg_D
lev = dimensions{1}.lev; % TODO : generalize to lev_D

%%
% Catch for TODO

% if nDims>1
%     for d=2:nDims
%         assert(dimensions{d}.lev==lev);
%         assert(dimensions{d}.deg==deg);
%     end
% end

%%
% TODO : surely this size depends on the parameters in the matrix_plot_D
% routine? i.e., the grid upon which we decide to evaluate the solution?

for d=1:nDims
%     dof_1D_FG(d) = deg*2^(dimensions{d}.lev);
    dof_1D_FG(d) = numel(Meval_D(:,1));
end
num_pts = prod(dof_1D_FG);

f_rSpace = sparse(num_pts,1);

for i=1:num_elements
    
    if opts.use_oldhash
        ll=hash_table{i};
    else
    end
    
    clear kronMatList;
    for d=1:nDims
        
        if opts.use_oldhash
            ID = ll(nDims*2+d);
        else
            IDlev = hash_table.elements.lev_p1(hash_table.elements_idx(i),d)-1;
            IDpos = hash_table.elements.pos_p1(hash_table.elements_idx(i),d)-1;
            ID = lev_cell_to_1D_index(IDlev,IDpos);
        end
        index_D = [(ID-1)*deg+1 : ID*deg];
        
        thisMat = Meval_D{d};
        thisMat1 = thisMat(:,index_D);
        
        kron_mat_list{d} = thisMat1; % Build kron list
    end
   
    element_ii = deg^nDims*(i-1)+1:deg^nDims*i;
    
    X = f_wSpace(element_ii);

    if use_kronmultd
        Y = kron_multd(nDims,kron_mat_list,X);
    else
        Y = kron_multd_full(nDims,kron_mat_list,X);
    end
        
    f_rSpace = f_rSpace + Y;
    
end

end



