function [V,H] = myarnoldi(A,v,M)
%given matrix A, finds a basis V for the M-dimensional Krylov subspace
%generated by {v,Av,A^2v,...,A^(M-1)v}
%H is the upper Hessenberg matrix.

[D,~] = size(A);
H = zeros(M+1,M);
V = zeros(D,M+1);

V(:,1) = v/norm(v);

for k=2:M+1
    V(:,k) = A*V(:,k-1);
    for j=1:k-1
        H(j,k-1) = V(:,j)'*V(:,k);
        V(:,k)=V(:,k)-H(j,k-1)*V(:,j); %orthogonalization
    end
    H(k,k-1)=norm(V(:,k));
    V(:,k) = V(:,k)/H(k,k-1);
end

V = V([1:D],[1:M]);
H = H([1:M],[1:M]);

end