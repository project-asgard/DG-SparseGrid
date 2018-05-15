function Y = kronmult1(A1, X )
% Y = kronmult1(A1, X )
nrow1 = size(A1,1);
ncol1 = size(A1,2);


Y = A1 * reshape(X, ncol1, numel(X)/ncol1 );

