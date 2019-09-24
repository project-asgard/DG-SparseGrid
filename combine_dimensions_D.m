function [fval] = combine_dimensions_D (deg, function_D, time_multiplier, hash_table, use_oldhash)

% Combine (via kron product) a set of 1D multiwavelet transforms to form
% the higher D sparse-grid multiwavelet representation.

num_dims = numel(function_D);

if use_oldhash
    num_elements = numel(hash_table);
else
    num_elements = numel(hash_table.elements_idx);
end

fval = zeros(deg^num_dims * num_elements,1);

for i=1:num_elements
    
    %%
    % Kron product approach
    
    for d=1:num_dims
        
        if use_oldhash
            ll=hash_table{i};
            idx_1D = ll(num_dims*2+d);
        else
            lev = hash_table.elements.lev_p1(hash_table.elements_idx(i),d)-1;
            pos = hash_table.elements.pos_p1(hash_table.elements_idx(i),d)-1;
            idx_1D = lev_cell_to_1D_index(lev,pos);
        end
        
        index_D = [(idx_1D-1)*deg+1 : idx_1D*deg];
        kron_mat_list{d} = function_D{d}(index_D);
    end
    
    B = time_multiplier;
    
    use_krond = 1;
    if use_krond
        A = krond (num_dims, kron_mat_list);
    else
        A = 1;
        for d=1:num_dims
            A = kron (A, kron_mat_list{d});
        end
    end
    
    tmp = A * B;
    
    Index = deg^num_dims*(i-1)+1:deg^num_dims*i;
    fval(Index,1) = fval(Index,1) + tmp(:);
    
end

end


