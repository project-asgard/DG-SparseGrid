function [f_nD] = singleD_to_multiD(num_dims,f_realspace,nodes)
if num_dims == 1
    n1 = numel(nodes{1});
    f_nD = f_realspace;
elseif num_dims == 2
    n1 = numel(nodes{1});
    n2 = numel(nodes{2});
    f_nD = reshape(f_realspace,n2,n1);
elseif num_dims == 3
    n1 = numel(nodes{1});
    n2 = numel(nodes{2});
    n3 = numel(nodes{3});
    f_nD = reshape(f_realspace,n3,n2,n1);
else
    error('Save output for num_dimensions >3 not yet implemented');
end

% % Remove duplicates
% 
% if strcmp(opts.output_grid,'fixed')
%     f_nD_nodups = ...
%         remove_duplicates(num_dims,f_nD,nodes_nodups,nodes_count);
%     f_nD = f_nD_nodups;
% end

end