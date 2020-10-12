function f_realspace_nD_fixed_grid_nodups = ...
    remove_duplicates(num_dimensions,f_realspace_nD_fixed_grid,nodes_fixed_grid_nodups,nodes_count)

if num_dimensions == 1
    dim = 1;
    ii = 1;
    for i=1:numel(nodes_fixed_grid_nodups{dim})
        if nodes_count{dim}(i) == 1
            f_realspace_nD_fixed_grid_nodups(i) = ...
                f_realspace_nD_fixed_grid(ii);
            ii = ii + 1;
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
    assert(sum(nodes_count{dim})==numel(f_realspace_nD_fixed_grid(1,:)));
    for i=1:numel(nodes_fixed_grid_nodups{dim})
        num_dups = nodes_count{dim}(i);
        if num_dups == 1
            f_realspace_nD_fixed_grid_nodups(:,i) = ...
                f_realspace_nD_fixed_grid(:,ii);
            ii = ii + 1;
        else
            f_realspace_nD_fixed_grid_nodups(:,i) = ...
                sum(f_realspace_nD_fixed_grid(:,ii:ii+num_dups-1),2)/num_dups;
            ii = ii + num_dups;
        end
    end
    % compress / remove duplicates from dim 2
    f_realspace_nD_fixed_grid = f_realspace_nD_fixed_grid_nodups;
    clear f_realspace_nD_fixed_grid_nodups;
    dim = 2;
    ii = 1;
    assert(sum(nodes_count{dim})==numel(f_realspace_nD_fixed_grid(:,1)));
    for i=1:numel(nodes_fixed_grid_nodups{dim})
        num_dups = nodes_count{dim}(i);
        if num_dups == 1
%             disp(num2str([i,ii,num_dups]));
            
            f_realspace_nD_fixed_grid_nodups(i,:) = ...
                f_realspace_nD_fixed_grid(ii,:);
            ii = ii + 1;
        else
%             disp(num2str([i,ii,num_dups]));
            
            f_realspace_nD_fixed_grid_nodups(i,:) = ...
                sum(f_realspace_nD_fixed_grid(ii:ii+num_dups-1,:),1)/num_dups;
            ii = ii + num_dups;
        end
    end
else    
    error('ERROR: remove_duplicates not implemented for dimensions > 2');
end

end