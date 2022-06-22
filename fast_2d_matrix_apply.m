function [f1] = fast_2d_matrix_apply(opts,pde,A_data,f0,imex_flag)
%Applys stiffness matrix-vector multiply f1=Af0 using the identity
%
% (B\kron A)*vec(X) = vec(AXB') where vec(X) = X(:);
%
persistent perm iperm pvec A_data_size

%Preallocate these
% if isempty(perm)
%     [perm,iperm,pvec] = sg_to_fg_mapping_2d(pde,opts,A_data);
% end
% if opts.adapt
%     [perm,iperm,pvec] = sg_to_fg_mapping_2d(pde,opts,A_data);
% end
if isempty(perm)
    [perm,iperm,pvec] = sg_to_fg_mapping_2d(pde,opts,A_data);
    A_data_size = numel(A_data.element_global_row_index);
end

if A_data_size ~= numel(A_data.element_global_row_index) %recompute
    [perm,iperm,pvec] = sg_to_fg_mapping_2d(pde,opts,A_data);
    A_data_size = numel(A_data.element_global_row_index);
end

%
if strcmp(opts.timestep_method,'IMEX')
    assert(strcmp(imex_flag,'E') || strcmp(imex_flag,'I'))
    term_idx = [];
    for i=1:numel(pde.terms)
        if strcmp(pde.terms{i}.imex,imex_flag)
            term_idx = [term_idx i];
        end
    end
else
    term_idx = 1:numel(pde.terms);
end

%Get dimensions
%n = size(pde.terms{1}.terms_1D{1}.mat,1);
%m = size(pde.terms{1}.terms_1D{2}.mat,1);
n = 2^pde.dimensions{1}.lev*opts.deg;
m = 2^pde.dimensions{2}.lev*opts.deg;

%Convert to standard FG index
f0_F = zeros(numel(perm),1);
f0_F(pvec) = f0(perm(pvec));

%Convert to matrix
f0_M = reshape(f0_F,m,n);

%Evaluate sum of krons
f1_M = zeros(size(f0_M));
for i=1:numel(term_idx)
    f1_M = f1_M + pde.terms{term_idx(i)}.terms_1D{2}.mat(1:m,1:m)*f0_M*pde.terms{term_idx(i)}.terms_1D{1}.mat(1:n,1:n)';
end

%Convert to vector
f1_F = f1_M(:);

%Back to SG-like index
f1 = f1_F(iperm);

end

