function [f_nD] = singleD_to_multiD(num_dims,fval_realspace,nodes)
if num_dims == 1
    n1 = numel(nodes{1});
    f_nD = fval_realspace;
elseif num_dims == 2
    n1 = numel(nodes{1});
    n2 = numel(nodes{2});
    f_nD = reshape(fval_realspace,n2,n1);
elseif num_dims == 3
    n1 = numel(nodes{1});
    n2 = numel(nodes{2});
    n3 = numel(nodes{3});
    f_nD = reshape(fval_realspace,n3,n2,n1);
else
    error('Save output for num_dimensions >3 not yet implemented');
end
end