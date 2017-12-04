function Y = kron_mult5(A1,A2,A3,A4,A5,X)
% Y = kron_mult3(A1,A2,A3,A4,A5,X)
%
Acell{1} = A1;
Acell{2} = A2;
Acell{3} = A3;
Acell{4} = A4;
Acell{5} = A5;
nkron = 5;
Y = kron_multd(nkron,Acell,X);

