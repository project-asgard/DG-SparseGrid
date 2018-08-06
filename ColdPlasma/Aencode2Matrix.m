function Mat = Aencode2Matrix(A_encode,dof)

Mat = sparse(dof,dof);

for i=1:size(A_encode,2)
    tmpA=A_encode{i}.A1;
    tmpB=A_encode{i}.A2;
    tmpC=A_encode{i}.A3;
    
    
    IndexI=A_encode{i}.IndexI;
    IndexJ=A_encode{i}.IndexJ;
    nI = size(IndexI,1);
    nJ = size(IndexJ,1);
    
    
    val = kron(tmpA,kron(tmpB,tmpC));
    
    
%     Mat = Mat+sparse(ones(nJ,1)*IndexI',IndexJ*ones(1,nI),val,dof,dof);
%     spy(Mat)
    Mat = Mat+sparse(IndexI*ones(1,nJ),ones(nI,1)*IndexJ',val,dof,dof);
end