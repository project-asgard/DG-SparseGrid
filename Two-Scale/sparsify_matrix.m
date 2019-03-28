function A_sparse = sparsify_matrix(A)
% A_sparse = sparsify_matrix(A)
%
% generate a sparse matrix representation
%
if (ismatrix(A) && ~issparse(A)),
        [ii,jj,aa] = find(A);
        A_sparse = sparse(ii,jj,aa, size(A,1), size(A,2));
else
     A_sparse = A;
end;

end
    
