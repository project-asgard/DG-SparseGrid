function [V,H] = arnoldi(A,q1,M)
%ARNOLDI    Arnoldi iteration
%   [Q,H] = ARNOLDI(A,q1,M) carries out M iterations of the
%   Arnoldi iteration with N-by-N matrix A and starting vector q1
%   (which need not have unit 2-norm).  For M < N it produces
%   an N-by-(M+1) matrix Q with orthonormal columns and an
%   (M+1)-by-M upper Hessenberg matrix H such that
%   A*Q(:,1:M) = Q(:,1:M)*H(1:M,1:M) + H(M+1,M)*Q(:,M+1)*E_M',
%   where E_M is the M'th column of the M-by-M identity matrix.

n = length(A);
if nargin < 3, M = n; end
q1 = q1/norm(q1);
V = zeros(n,M); V(:,1) = q1;
H = zeros(min(M+1,M),n);

for k=1:M
    z = A*V(:,k);
    for i=1:k
        H(i,k) = V(:,i)'*z;
        z = z - H(i,k)*V(:,i);
    end
    if k < n
       H(k+1,k) = norm(z);
       if H(k+1,k) == 0, return, end
       V(:,k+1) = z/H(k+1,k);
   end
end