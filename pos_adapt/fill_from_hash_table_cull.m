function [f,hash_table,A_data] = fill_from_hash_table_cull(pde,opts,hash_table_old,hash_table,f_old,f_bc,elements_removed,A_data,data)
%After refinement, the function needs to be backfilled. 
%This algorithm does that.  

num_ele_old = numel(hash_table_old.elements_idx);
num_ele = numel(hash_table.elements_idx);

num_dims = numel(pde.dimensions);
deg = opts.deg;

cell_dof = deg^num_dims;

%Assert that the new hash_table is just a refinement of the old one
assert(sum(hash_table_old.elements_idx-hash_table.elements_idx(1:num_ele_old)) == 0) %This is an integer sum so setting to zero is okay.
assert(num_ele_old*cell_dof == numel(f_old));


if nargin < 6 %Element data added need to be zero

    f = zeros(num_ele*cell_dof,1);
    f(1:num_ele_old*cell_dof) = f_old;

else %We only need to add element data from f_bc
    
    start = num_ele_old*cell_dof;
    
    cells_to_add = hash_table.elements_idx(num_ele_old+1:num_ele);
    cells_not_added = [];
    %Determine which elements from new hash we removed from f in the
    %adaptive culling
    [I,loc] = ismember(cells_to_add,elements_removed(2,:));
    
    f = zeros(start+(sum(I)*cell_dof),1);
    f(1:num_ele_old*cell_dof) = f_old;
    
    count = 0;
    for i=1:numel(cells_to_add)
        if I(i)
            count = count + 1;
            %Get idx of cell
            idx = elements_removed(1,loc(i));
            f(start+(count-1)*cell_dof+1:start+count*cell_dof) = f_bc((idx-1)*cell_dof+1:idx*cell_dof);
            %f(start+(count-1)*cell_dof+1) = f_bc((idx-1)*cell_dof+1);
        else
            %Remove element from hash_table
            cells_not_added = [cells_not_added,num_ele_old+i];
            hash_table.elements.lev_p1(hash_table.elements_idx(num_ele_old+i),:) = 0;
            hash_table.elements.pos_p1(hash_table.elements_idx(num_ele_old+i),:) = 0;
            hash_table.elements.type(hash_table.elements_idx(num_ele_old+i))= 0;            
        end
    end
    
    %Remove hash_table elements that didnt get added
    if ~isempty(cells_not_added)
         for ele=cells_not_added
             hash_table.elements.lev_p1(hash_table.elements_idx(ele),:) = 0;
             hash_table.elements.pos_p1(hash_table.elements_idx(ele),:) = 0;
             hash_table.elements.type(hash_table.elements_idx(ele))= 0;
        end
        hash_table.elements_idx(cells_not_added) = []; 
        A_data = global_matrix(pde,opts,hash_table);
    end
    
    [perm,~,pvec] = sg_to_fg_mapping_2d(pde,opts,A_data);
    f_FG = zeros(size(pvec));
    f_FG(pvec) = f(perm(pvec));
    Q_new = data.B*f_FG;
end


end

