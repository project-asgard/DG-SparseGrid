function Y = kron_mult6(A1,A2,A3,A4,A5,A6,X)
% Y = kron_mult6(A1,A2,A3,A4,A5,A6,X)
%
Acell{1} = A1;
Acell{2} = A2;
Acell{3} = A3;
Acell{4} = A4;
Acell{5} = A5;
Acell{6} = A6;
nkron = 6;
Y = kron_multd(nkron,Acell,X);

