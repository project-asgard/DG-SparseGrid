function flops = kron_cost_fixed( rc )
%
% flops = kron_cost_fixed( rc )
%
% return flops for fixed strategy
% (An * X) * {kron(A1...Anm1)}^t
% or
% X is n2 by n1
% Z = An * X   (m2 by n1)
% Zt = transpose(Z) is n1 by m2
% kron(A1..Anm1)*Zt
% An is m2 by n2,  kron(A1...Anm1) is m1 by n1
% --------------------------------------
nkron = size(rc,2);
n = nkron;


if (nkron == 1),
  % --------------------
  % just a single matrix
  % --------------------
  nrowA = rc(1,1);
  ncolA = rc(2,1);
  flops = 2.0*nrowA * ncolA;
else
 % -------------
 % use recursion
 % -------------
 nm1 = n - 1;
 m1 = prod( rc(1,1:nm1) );
 n1 = prod( rc(2,1:nm1) );
 m2 = rc(1,n);
 n2 = rc(1,n);

 flops_Z = 2.0*m2 * n2 * n1;
 flops1 = kron_cost_fixed( rc(1:2,1:(nm1)) );
 flops = flops1 * m2 + flops_Z;
end;
return;
end;

   

