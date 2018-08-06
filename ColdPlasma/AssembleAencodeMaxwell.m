function A = AssembleAencodeMaxwell(IndexI,IndexJ,GG,G,J,H,K,L,Q,w2,dof,row1,row2,row3,col1,col2,col3,alpha)
% Assemble the global matrix
% dof: dof for any of E1(x,y,z),E2(x,y,z),E3(x,y,z)

GT=G';

nx1 = length(row1);ny1 = length(col1);
nx2 = length(row2);ny2 = length(col2);
nx3 = length(row3);ny3 = length(col3);

dof_1D=size(G,1);

I1 = speye(dof_1D,dof_1D);
I2 = speye(dof_1D,dof_1D);
I3 = speye(dof_1D,dof_1D);
% A11
count = 1;
A{count}.IndexI = IndexI';
A{count}.IndexJ = IndexJ';

A{count}.A1=I1(row1,col1);%speye(nx1,ny1);
A{count}.A2=GG(row2,col2)-J(row2,col2)-L(row2,col2)+alpha*Q(row2,col2)-w2*I2(row2,col2);%
A{count}.A3=I3(row3,col3);%
%     full(sparse(IndexI',IndexJ',kron(II,kron(A{count}.A2,II))))

count = count+1;
A{count}.IndexI = IndexI';
A{count}.IndexJ = IndexJ';

A{count}.A1=I1(row1,col1);%speye(nx1,ny1);
A{count}.A2=I2(row2,col2);%
A{count}.A3=GG(row3,col3)-J(row3,col3)-L(row3,col3)+alpha*Q(row3,col3);

%     full(sparse(IndexI',IndexJ',kron(II,kron(II,A{count}.A3))))

% A12
count = count+1;
A{count}.IndexI = IndexI';
A{count}.IndexJ = dof+IndexJ';

A{count}.A1=-G(row1,col1);
A{count}.A2= GT(row2,col2);%G(row2,col2)';
A{count}.A3=I3(row3,col3);%

count = count+1;
A{count}.IndexI = IndexI';
A{count}.IndexJ = dof+IndexJ';

A{count}.A1=H(row1,col1);
A{count}.A2=GT(row2,col2);%G(row2,col2)';
A{count}.A3=I3(row3,col3);%

count = count+1;
A{count}.IndexI = IndexI';
A{count}.IndexJ = dof+IndexJ';

A{count}.A1=G(row1,col1);
A{count}.A2=K(row2,col2);%H';
A{count}.A3=I3(row3,col3);%

% A13
count = count+1;
A{count}.IndexI = IndexI';
A{count}.IndexJ = 2*dof+IndexJ';

A{count}.A1=-G(row1,col1);
A{count}.A2=I2(row2,col2);%
A{count}.A3=GT(row3,col3);%G(row3,col3)';

count = count+1;
A{count}.IndexI = IndexI';
A{count}.IndexJ = 2*dof+IndexJ';

A{count}.A1=H(row1,col1);
A{count}.A2=I2(row2,col2);%
A{count}.A3=GT(row3,col3);%G(row3,col3)';

count = count+1;
A{count}.IndexI = IndexI';
A{count}.IndexJ = 2*dof+IndexJ';

A{count}.A1=G(row1,col1);
A{count}.A2=I2(row2,col2);%
A{count}.A3=K(row3,col3);%H';

% A21
count = count+1;
A{count}.IndexI = dof+IndexI';
A{count}.IndexJ = IndexJ';

A{count}.A1=-GT(row1,col1);%G(row1,col1)';
A{count}.A2=G(row2,col2);
A{count}.A3=I3(row3,col3);%

count = count+1;
A{count}.IndexI = dof+IndexI';
A{count}.IndexJ = IndexJ';

A{count}.A1=GT(row1,col1);%G(row1,col1)';
A{count}.A2=H(row2,col2);
A{count}.A3=I3(row3,col3);%

count = count+1;
A{count}.IndexI = dof+IndexI';
A{count}.IndexJ = IndexJ';

A{count}.A1=K(row1,col1);%H';
A{count}.A2=G(row2,col2);
A{count}.A3=I3(row3,col3);%

% A22
count = count+1;
A{count}.IndexI = dof+IndexI';
A{count}.IndexJ = dof+IndexJ';

A{count}.A1=GG(row1,col1)-J(row1,col1)-L(row1,col1)+alpha*Q(row1,col1)-w2*I1(row1,col1);%speye(nx1,ny1);
A{count}.A2=I2(row2,col2);%
A{count}.A3=I3(row3,col3);%

count = count+1;
A{count}.IndexI = dof+IndexI';
A{count}.IndexJ = dof+IndexJ';

A{count}.A1=I1(row1,col1);%
A{count}.A2=I2(row2,col2);%
A{count}.A3=GG(row3,col3)-J(row3,col3)-L(row3,col3)+alpha*Q(row3,col3);

% A23
count = count+1;
A{count}.IndexI = dof+IndexI';
A{count}.IndexJ = 2*dof+IndexJ';

A{count}.A1=-I1(row1,col1);%
A{count}.A2=G(row2,col2);
A{count}.A3=GT(row3,col3);%G(row3,col3)';

count = count+1;
A{count}.IndexI = dof+IndexI';
A{count}.IndexJ = 2*dof+IndexJ';

A{count}.A1=I1(row1,col1);%
A{count}.A2=H(row2,col2);
A{count}.A3=GT(row3,col3);%G(row3,col3)';

count = count+1;
A{count}.IndexI = dof+IndexI';
A{count}.IndexJ = 2*dof+IndexJ';

A{count}.A1=I1(row1,col1);%
A{count}.A2=G(row2,col2);
A{count}.A3=K(row3,col3);%H';

% A31
count = count+1;
A{count}.IndexI = 2*dof+IndexI';
A{count}.IndexJ = IndexJ';

A{count}.A1=-GT(row1,col1);%G(row1,col1)';
A{count}.A2=I2(row2,col2);%
A{count}.A3=G(row3,col3);

count = count+1;
A{count}.IndexI = 2*dof+IndexI';
A{count}.IndexJ = IndexJ';

A{count}.A1=GT(row1,col1);%G(row1,col1)';
A{count}.A2=I2(row2,col2);%
A{count}.A3=H(row3,col3);

count = count+1;
A{count}.IndexI = 2*dof+IndexI';
A{count}.IndexJ = IndexJ';

A{count}.A1=K(row1,col1);%H';
A{count}.A2=I2(row2,col2);%
A{count}.A3=G(row3,col3);

% A32
count = count+1;
A{count}.IndexI = 2*dof+IndexI';
A{count}.IndexJ = dof+IndexJ';

A{count}.A1=-I1(row1,col1);%
A{count}.A2=GT(row2,col2);%G(row2,col2)';
A{count}.A3=G(row3,col3);

count = count+1;
A{count}.IndexI = 2*dof+IndexI';
A{count}.IndexJ = dof+IndexJ';

A{count}.A1=I1(row1,col1);%
A{count}.A2=GT(row2,col2);%G(row2,col2)';
A{count}.A3=H(row3,col3);

count = count+1;
A{count}.IndexI = 2*dof+IndexI';
A{count}.IndexJ = dof+IndexJ';

A{count}.A1=I1(row1,col1);%
A{count}.A2=K(row2,col2);%H';
A{count}.A3=G(row3,col3);

% A33
count = count+1;
A{count}.IndexI = 2*dof+IndexI';
A{count}.IndexJ = 2*dof+IndexJ';

A{count}.A1=GG(row1,col1)-J(row1,col1)-L(row1,col1)+alpha*Q(row1,col1)-w2*I1(row1,col1);%speye(nx1,ny1);
A{count}.A2=I2(row2,col2);%
A{count}.A3=I3(row3,col3);%

count = count+1;
A{count}.IndexI = 2*dof+IndexI';
A{count}.IndexJ = 2*dof+IndexJ';

A{count}.A1=I1(row1,col1);%
A{count}.A2=GG(row2,col2)-J(row2,col2)-L(row2,col2)+alpha*Q(row2,col2);
A{count}.A3=I3(row3,col3);%