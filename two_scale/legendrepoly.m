function F = legendrepoly( ndeg, x, do_normalize_in )
% function F = legendrepoly( ndeg, x )
% evaluate orthogonal legendre polynomial of degree 0, 1, ..., ndeg
% at positions x(1), x(2), ... x(n)
% F(i,j) = P_{j-1}( x(i) ),   
% where P_j(x) is legendre polynomial of degree j
%
do_normalize = 0;
if (nargin >= 3),
   do_normalize = do_normalize_in;
end;

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

if (do_normalize),
% ----------------------------------------------------
% integral( P_n(x) * P_n(x), over x=-1..1 ) = 2/(2*n +1)
% thus rescale Phat_n(x) =  P_n(x)/sqrt( 2/(2*n+1) )
% integral( Phat_n(x) * Phat_n(x), over x=-1..1 ) = 1
% ----------------------------------------------------
for L=0:ndeg,
   dnorm = 2/(2*L+1);
   dscale = sqrt( dnorm );
   F(1:n, ioff+L) = F(1:n,ioff+L)/dscale;
end;

end

