function Y = kron_mult3(A1,A2,A3,X)
% Y = kron_mult3(A1,A2,A3,X)
%
Acell{1} = A1;
Acell{2} = A2;
Acell{3} = A3;
nkron = 3;
Y = kron_multd(nkron,Acell,X);

