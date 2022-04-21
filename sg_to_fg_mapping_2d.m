function [perm,iperm,pvec] = sg_to_fg_mapping_2d(md_term,opts,A_data)
% Maps the sparse_grid index to standard tensor product full_grid index.
% This to allow multiplication by the full-grid kronecker sum stiffness
% matrix.  

% perm  -> SG-type index to FG-standard index
% pvec  -> logical vector showing which elements are lit.  Used for SG->FG
%            conversion.
% iperm -> inverse mapping

% This reindexing allows us to use the standard sum of kroncker products 
%  stiffness matrix in lieu of the construction used in apply_A.  Using the
%  conversion to full-grid and back is also much faster.

% Works for adaptiviity

% Only works in 2D for now, but can be easily changed for higher diminsions

% Function assumes A_data.element_global_row_index is 1:sparse_dof.
%  --Does not check for this

% USAGE:  If the matrix A and vector f are in the SG-index, then
%   we can compute A*f=f1 using A_F by the sequence
%     f_F = zeros(size(A_F,1),1);
%     f_F(pvec) = f(perm(pvec));
%     Af_F = A_F*f_F;
%     f1  = Af_F(iperm);
% where A_F is the standard full-grid sum of kronecker products.
% Specifically, A_F is generated by the code:
%   A_F = 0;
%   for i=1:numel(pde.terms)
%       A_F = A_F + kron(pde.terms{i}.terms_1D{1}.mat,...
%                        pde.terms{i}.terms_1D{2}.mat);
%   end
%
% The SG-index stiffness matrix and be computed from the full grid matrix
% by using:
%         A = A_F(iperm,iperm);

%Checking for only 2 dimensions
assert(numel(md_term.terms_1D) == 2);

deg = opts.deg;

num_sparse_ele = numel(A_data.element_local_index_D{1}); %Number of sparse elements
iperm = zeros(num_sparse_ele*deg^2,1);

%Get FG_dof for each direction
max_x_dof = size(md_term.terms_1D{1}.mat,1);
max_y_dof = size(md_term.terms_1D{2}.mat,1);
perm = nan(max_x_dof*max_y_dof,1);

idx = 1; %Keep track of index
for i=1:num_sparse_ele
    %First compute the x dim FG dof
    dof_x = (A_data.element_local_index_D{1}(i)-1)*deg + (1:deg);
    %Now y dim FG dof
    dof_y = (A_data.element_local_index_D{2}(i)-1)*deg + (1:deg);
    
    %Mapping is (i,j) -> (i-1)*max_y_dof+j;
    
    for deg_x=1:deg
        for deg_y=1:deg
            iperm(idx) = (dof_x(deg_x)-1)*max_y_dof+dof_y(deg_y);
            idx = idx + 1;
        end
    end
end


perm(iperm) = 1:numel(iperm);
pvec = ~isnan(perm);


end
