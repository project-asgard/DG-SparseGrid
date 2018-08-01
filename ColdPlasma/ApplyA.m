function ftmp = ApplyA(b)
% A*b
global A_encode

dof = size(b,1);
ftmp=sparse(dof,1);

nkron = 3;

for i=1:size(A_encode,2)
tmpA=A_encode{i}.A1;
tmpB=A_encode{i}.A2;
tmpC=A_encode{i}.A3;
Acell{1} = tmpA;
Acell{2} = tmpB;
Acell{3} = tmpC;

IndexI=A_encode{i}.IndexI;
IndexJ=A_encode{i}.IndexJ;
ftmp(IndexI)=ftmp(IndexI)+kron_multd(nkron,Acell,b(IndexJ));
clear Acell
end
