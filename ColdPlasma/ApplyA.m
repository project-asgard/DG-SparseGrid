function ftmp = ApplyA(f)
% A*b
global A_encode 

dof = size(f,1);
ftmp=zeros(dof,1);


for i=1:size(A_encode,2)
tmpA=A_encode{i}.A1;
tmpB=A_encode{i}.A2;
tmpC=A_encode{i}.A3;


IndexI=A_encode{i}.IndexI;
IndexJ=A_encode{i}.IndexJ;


val = kron_mult3(tmpA,tmpB,tmpC,f(IndexJ));
ftmp(IndexI)=ftmp(IndexI)+val;

end


