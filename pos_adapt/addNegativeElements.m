function [hash_table,A_data] = addNegativeElements(pde,opts,hash_table,Q,pos_tol)
%Q represents the constant polynomial coefficients of the function on the
%full grid.  If the value is negative, this function adds the hierarchical
%basis functions to the new hash to hopefully perserve positivity of the
%cell averages.

persistent Ix Iv


assert(numel(pde.dimensions) == 2);
num_dims = numel(pde.dimensions);

n_v = uint64(2^pde.dimensions{2}.lev);
lev_x = pde.dimensions{1}.lev;
lev_v = pde.dimensions{2}.lev;
max_lev = max([lev_x,lev_v]);

if isempty(Ix)
    [Ix,~] = find(OperatorTwoScale_wavelet2(1,lev_x));
    [Iv,~] = find(OperatorTwoScale_wavelet2(1,lev_v));
end

[I] = uint64(find(Q < pos_tol));

x_idx = idivide(I-1,n_v) + 1;
v_idx = mod(I-1,n_v) + 1;

pre_elements_to_add = zeros(numel(x_idx)*(lev_x+1)*(lev_v+1),2);
%Need to go from realspace to hierarchical space
for i=1:numel(x_idx)
    x_hier = Ix((x_idx(i)-1)*(lev_x+1)+1:x_idx(i)*(lev_x+1));
    v_hier = Iv((v_idx(i)-1)*(lev_v+1)+1:v_idx(i)*(lev_v+1));
    
    %Matlab magic to do list all combinations of these two vectors from
    % -------
    % https://www.mathworks.com/matlabcentral/answers/98191-how-can-i-obtain-all-possible-combinations-of-given-vectors-in-matlab
    % -------
    [A,B] = meshgrid(x_hier,v_hier);
    C = cat(2,A',B');
    pre_elements_to_add( (i-1)*(lev_x+1)*(lev_v+1)+1:i*(lev_x+1)*(lev_v+1),:) = reshape(C,[],2);
end

%Get global md index of each element
idx_vec = pre_elements_to_add(:,1) + (pre_elements_to_add(:,2)-1)*2^max_lev;

%Get lev and position of each element
hash_old = hash_table;

num_elements = numel(hash_table.elements_idx);
num_elements_added = 0;
for i=1:numel(idx_vec)
    idx = idx_vec(i);
    
    if hash_table.elements.type(idx) == 0 % element not already enabled and level does not grow
            
        num_elements_added = num_elements_added + 1;
        position_in_elements_idx = num_elements+num_elements_added;
        hash_table.elements_idx(position_in_elements_idx) = idx; % Extend element list
        
        [lev_vec, pos_vec] = md_idx_to_lev_pos(num_dims, opts.max_lev, idx);

        hash_table.elements.lev_p1(idx,:) = lev_vec+1; % NOTE : have to start lev  index from 1 for sparse storage
        hash_table.elements.pos_p1(idx,:) = pos_vec+1; % NOTE : have to start cell index from 1 for sparse storage
        hash_table.elements.type(idx) = 1;

    end
    
end

A_data = global_matrix(pde,opts,hash_table);

end

