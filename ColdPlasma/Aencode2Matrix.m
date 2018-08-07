function Mat = Aencode2Matrix(A_encode,dof)

Mat = sparse(dof,dof);

for i=1:size(A_encode,2)
    tmpA=A_encode{i}.A1;
    tmpB=A_encode{i}.A2;
    tmpC=A_encode{i}.A3;
    
    
    IndexI=A_encode{i}.IndexI;
    IndexJ=A_encode{i}.IndexJ;
    nI = length(IndexI);
    nJ = length(IndexJ);
    
    
    val = kron(kron(tmpA,tmpB),tmpC);
%     ones(nJ,1)*IndexI
%     
%     IndexJ'*ones(1,nI)
%     Mat = Mat + sparse(IndexI'*ones(1,nJ),ones(nI,1)*IndexJ,val,dof,dof);
%     spy(Mat)
    Mat = Mat+sparse(IndexI*ones(1,nJ),ones(nI,1)*IndexJ',val,dof,dof);
end