function ftmp = ApplyA_imag(f)
% A*b
global A_encode 

dof = size(f,1);
ftmp=zeros(dof,1);


for i=1:size(A_encode,2)
tmpA=A_encode{i}.A;
tmpB=A_encode{i}.B;


IndexI=A_encode{i}.IndexI;
IndexJ=A_encode{i}.IndexJ;

val = kronmult2(imag(tmpA),imag(tmpB),f(IndexJ));
ftmp(IndexI)=ftmp(IndexI)+val;

end


