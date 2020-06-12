function f_realspace_nD_fixed_grid_nodups = ...
    remove_duplicates(num_dimensions,f_realspace_nD_fixed_grid,nodes_fixed_grid_nodups,nodes_count)

if num_dimensions == 1
    dim = 1;
    ii = 1;
    for i=1:numel(nodes_fixed_grid_nodups{dim})
        if nodes_count{dim} == 1
            f_realspace_nD_fixed_grid_nodups(i) = ...
                f_realspace_nD_fixed_grid(ii);
            ii = ii + i;
        else
            num_dups = nodes_count{1}(i);
            f_realspace_nD_fixed_grid_nodups(i) = ...
                sum(f_realspace_nD_fixed_grid(ii:ii+num_dups-1),dim)/num_dups;
            ii = ii + num_dups;
        end
    end
elseif num_dimensions == 2
    % compress / remove duplicates from dim 1
    clear f_realspace_nD_fixed_grid_nodups;
    dim = 1;
    ii = 1;
    for i=1:numel(nodes_fixed_grid_nodups{dim})
        if nodes_count{dim} == 1
            f_realspace_nD_fixed_grid_nodups(i,:) = ...
                f_realspace_nD_fixed_grid(ii,:);
            ii = ii + i;
        else
            num_dups = nodes_count{dim}(i);
            f_realspace_nD_fixed_grid_nodups(i,:) = ...
                sum(f_realspace_nD_fixed_grid(ii:ii+num_dups-1,:),1)/num_dups;
            ii = ii + num_dups;
        end
    end
    % compress / remove duplicates from dim 2
    f_realspace_nD_fixed_grid = f_realspace_nD_fixed_grid_nodups;
    clear f_realspace_nD_fixed_grid_nodups;
    dim = 2;
    ii = 1;
    for i=1:numel(nodes_fixed_grid_nodups{dim})
        if nodes_count{dim} == 1
            f_realspace_nD_fixed_grid_nodups(:,i) = ...
                f_realspace_nD_fixed_grid(:,ii);
            ii = ii + i;
        else
            num_dups = nodes_count{dim}(i);
            f_realspace_nD_fixed_grid_nodups(:,i) = ...
                sum(f_realspace_nD_fixed_grid(:,ii:ii+num_dups-1),dim)/num_dups;
            ii = ii + num_dups;
        end
    end
end

end