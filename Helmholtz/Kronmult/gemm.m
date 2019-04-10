function [C] = gemm( transA, transB, m,n,k, alpha, A, B, beta, C )
%
% [C] = gemm( transA, transB, m,n,k, alpha, A, B, beta, C )
%
% simulate dgemm
%
opA = A;
if ((transA == 'T') || (transA == 't')),
  opA = transpose(A);
end;
if ((transA == 'C') || (transA == 'c')),
  opA = conj(transpose(A));
end;

opB = B;
if ((transB == 'T') || (transB == 't')),
  opB = transpose(B);
end;
if ((transB == 'C') || (transB == 'c')),
  opB = conj(transpose(B));
end;

if (beta == 0),
  C(1:m,1:n) = 0;
else
  C(1:m,1:n) = beta * C(1:m,1:n);
end;

C(1:m,1:n) = C(1:m,1:n) + alpha * opA(1:m, 1:k) * opB(1:k, 1:n);

