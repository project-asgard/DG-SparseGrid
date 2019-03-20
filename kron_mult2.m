function Y = kron_mult2(A1,A2,X)
% Y = kron_mult2(A1,A2,X)
%
use_kron_multd = 0;
if (use_kron_multd),
 Acell{1} = A1;
 Acell{2} = A2;
 nkron = 2;
 Y = kron_multd(nkron,Acell,X);
else
  Y = kronmult2(A1,A2,X);
end;
