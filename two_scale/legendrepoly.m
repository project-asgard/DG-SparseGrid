function F = legendrepoly( ndeg, x )
% function F = legendrepoly( ndeg, x )
% evaluate orthogonal legendre polynomial of degree 0, 1, ..., ndeg
% at positions x(1), x(2), ... x(n)
% F(i,j) = P_{j-1}( x(i) ),   
% where P_j(x) is legendre polynomial of degree j
%
n = numel(x);
x = reshape(x,n,1);
F = zeros( n,ndeg+1);
ioff = 1;
F(1:n,ioff+0) = 1;
if (ndeg >= 1),
   F(1:n,ioff+1) = x;
end;

% ----------------------------------------------------
% Legendre polynomials satisfy the recurrence relation
%  (L+1)*P_{L+1}(x) = (2*L+1)*x*P_{L}(x) - L*P_{L-1}(x)
% ----------------------------------------------------
for L=1:ndeg-1,
      PLm1 = F(1:n,ioff+ (L-1));
      PL = F(1:n,ioff+L);

      PLp1 = ((2*L+1)*x(1:n).*PL(1:n) - L*PLm1(1:n))/(L+1);
      F(1:n,ioff+ (L+1) )  = PLp1(1:n);
end;

end

