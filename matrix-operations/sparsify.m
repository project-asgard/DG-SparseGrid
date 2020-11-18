function As = sparsify( A, tol_in)
% As = sparsify( A, tol )
% generate sparse version of matrix A
% accept A(i,j) if  abs(A(i,j)) > tol
% -----------------------------------
tol = eps;
if (nargin >= 2),
    tol = tol_in;
end;


[i,j,aij] = find(A);
accept = (abs(aij) > tol);
nrow = size(A,1);
ncol = size(A,2);
As = sparse( i(accept), j(accept), aij(accept), nrow,ncol);
