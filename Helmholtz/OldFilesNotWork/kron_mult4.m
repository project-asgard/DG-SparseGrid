function Y = kron_mult4(A1,A2,A3,A4,X)
% Y = kron_mult4(A1,A2,A3,A4,X)
%
Acell{1} = A1;
Acell{2} = A2;
Acell{3} = A3;
Acell{4} = A4;
nkron = 4;
Y = kron_multd(nkron,Acell,X);

