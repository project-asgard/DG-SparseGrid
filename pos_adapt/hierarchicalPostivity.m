function [hash_new,A_new] = hierarchicalPostivity(pde,opts,f,hash_table,A_data,pos_tol)
%Hierarchical adaptivity for positivity

persistent FG2DG Ix

%Need to create hierarchical FG to realspace DG constant transfer matrix
%for each level
pde_lev_vec = [pde.dimensions{1}.lev,pde.dimensions{2}.lev];
deg = opts.deg;

hash_new = hash_table;
A_new = A_data;


%Compute persistent values for speed
if isempty(FG2DG)
    FG2DG = cell(max(pde_lev_vec),1);
    Ix = cell(max(pde_lev_vec),1);
    
    
    for l=1:max(pde_lev_vec) %Level 0 is always in the space.
        lev_vec = min([l l; pde_lev_vec]);
        
        %Build restriction matrices to lower level space
        R_x = speye(2^lev_vec(1)*deg,2^pde_lev_vec(1)*deg);
        R_v = speye(2^lev_vec(2)*deg,2^pde_lev_vec(2)*deg);
        R   = kron(R_x,R_v);
        FMWT_x = OperatorTwoScale_wavelet2(deg,lev_vec(1));
        FMWT_v = OperatorTwoScale_wavelet2(deg,lev_vec(2));
        FMWT_2D = kron(FMWT_x',FMWT_v');
        fg2dg    = kronrealspace2DtoDG(lev_vec,deg,1)*FMWT_2D;
        FG2DG{l} = fg2dg*R*sqrt(2^(lev_vec(1)-1)*2^(lev_vec(2)-1));
        
        %Lastly construst wavelet passthrough for representing local
        %elements in the hierarchical wavelet space
        [Ix{l},~] = find(OperatorTwoScale_wavelet2(1,l));
    end
    
end

[perm,~,pvec] = sg_to_fg_mapping_2d(pde,opts,A_data);
%I = speye(numel(iperm));
%I_FG = sparse(size(pvec,1),numel(iperm)); I_FG(pvec,:) = I(perm(pvec),:);

f_FG = zeros(size(pvec));
f_FG(pvec) = f(perm(pvec));

%Map to cell averages on coarse grids
for l=1:max(pde_lev_vec)
    lev_vec = min([l l; pde_lev_vec]);
    
    Q = FG2DG{l}*f_FG;
    fprintf('Level %d.  min(Q) = %e\n',l,min(Q));
    
    %Add elements if negative on this level
    if min(Q) < pos_tol
        hash_new = addHierNegativeElements(lev_vec,pde_lev_vec,opts,hash_table,Q,Ix{lev_vec(1)},Ix{lev_vec(2)},pos_tol);
        A_new = global_matrix(pde,opts,hash_new);
        break
    end
end




end

function hash_table = addHierNegativeElements(lev_vec,pde_lev_vec,opts,hash_table,Q,Ix,Iv,pos_tol)
%Q represents the constant polynomial coefficients of the function on the
%full grid.  If the value is negative, this function adds the hierarchical
%basis functions to the new hash to hopefully perserve positivity of the
%cell averages.

num_dims = numel(lev_vec);

n_v = uint64(2^lev_vec(2));
lev_x = lev_vec(1);
lev_v = lev_vec(2);
max_lev = max(pde_lev_vec);

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

end


