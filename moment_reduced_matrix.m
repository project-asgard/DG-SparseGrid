function [M] = moment_reduced_matrix(opts,pde,A_data,hash_table,moment_idx)
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CURRENTLY ONLY WORKS FOR 1x-1v system
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% This function produces a matrix M that is the transformation matrix that
% maps sparse-grid functions in (x,v) to DG functions only in x.  This
% transformation is given by 
%          f_x(x) = \int_v f(x,v)*g(v) dv
% where g(v) is a moment function with index moment_idx
% 
% In 1x-1v systems, the output function is full grid, but this will break
% down in multiple dimensions.   

deg = opts.deg;
num_ele = numel(A_data.element_local_index_D{1});
x_dim = 1; %Hardcoded for now.  Needs to change in multi-dim
v_dim = 2; %Hardcoded for now.  Needs to change in multi-dim
g_vec = pde.moments{moment_idx}.fList{v_dim};



%Get row-dim of M
rows = max(A_data.element_local_index_D{x_dim}*deg);

%Containers for sparse matrix
I = zeros(deg^2*num_ele,1); %Square is hardcoded.  Needs to change 
J = zeros(deg^2*num_ele,1);
S = zeros(deg^2*num_ele,1);

count = 1;
k = 1:deg;
for i=1:num_ele
    for j=1:deg
        temp = ((A_data.element_local_index_D{x_dim}(i)-1)*deg + j) + 0*k;
        I(count:count+deg-1) = temp';
        
        temp = (i-1)*deg^2 + (j-1)*deg + k;
        J(count:count+deg-1) = temp';
        
        temp = g_vec((A_data.element_local_index_D{v_dim}(i)-1)*deg + k);
        S(count:count+deg-1) = temp';
        
        count = count + deg;
    end
end

M = sparse(I,J,S,rows,deg^2*num_ele);


end

