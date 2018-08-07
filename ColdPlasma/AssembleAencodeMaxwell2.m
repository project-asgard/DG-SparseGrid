function A = AssembleAencodeMaxwell(IndexI,IndexJ,GG,G,J,H,K,L,Q,w2,dof,row1,row2,row3,col1,col2,col3,alpha)
% Assemble the global matrix
% dof: dof for any of E1(x,y,z),E2(x,y,z),E3(x,y,z)
A={};

GT=G';

nx1 = length(row1);ny1 = length(col1); z1 = zeros(nx1,ny1);
nx2 = length(row2);ny2 = length(col2); z2 = zeros(nx2,ny2);
nx3 = length(row3);ny3 = length(col3); z3 = zeros(nx3,ny3);

dof_1D=size(G,1);

I1 = speye(dof_1D,dof_1D);
I2 = speye(dof_1D,dof_1D);
I3 = speye(dof_1D,dof_1D);

count = 0;

% A11
A1 = I1(row1,col1);
A2 = GG(row2,col2)-J(row2,col2)-L(row2,col2)+alpha*Q(row2,col2)-w2*I2(row2,col2);%
A3 = I3(row3,col3);
% if norm(A1-z1)>1e-4 && norm(A2-z2)>1e-4 && norm(A3-z3)>1e-4
    count = count+1;
A{count}.IndexI = IndexI;
A{count}.IndexJ = IndexJ;

A{count}.A1 = A1;%I1(row1,col1);%speye(nx1,ny1);
A{count}.A2 = A2;%GG(row2,col2)-J(row2,col2)-L(row2,col2)+alpha*Q(row2,col2)-w2*I2(row2,col2);%
A{count}.A3 = A3;%I3(row3,col3);%
% end
%     full(sparse(IndexI,IndexJ,kron(II,kron(A{count}.A2,II))))

A1 = I1(row1,col1);
A2 = I2(row2,col2);%
A3 = GG(row3,col3)-J(row3,col3)-L(row3,col3)+alpha*Q(row3,col3);
% if norm(A1-z1)>1e-4 && norm(A2-z2)>1e-4 && norm(A3-z3)>1e-4
count = count+1;
A{count}.IndexI = IndexI;
A{count}.IndexJ = IndexJ;

A{count}.A1=A1;%I1(row1,col1);%speye(nx1,ny1);
A{count}.A2=A2;%I2(row2,col2);%
A{count}.A3=A3;%GG(row3,col3)-J(row3,col3)-L(row3,col3)+alpha*Q(row3,col3);
% end
%     full(sparse(IndexI,IndexJ,kron(II,kron(II,A{count}.A3))))

% A12
A1 = -G(row1,col1);
A2 = GT(row2,col2);%G(row2,col2)';
A3 = I3(row3,col3);%
% if norm(A1-z1)>1e-4 && norm(A2-z2)>1e-4 && norm(A3-z3)>1e-4
count = count+1;
A{count}.IndexI = IndexI;
A{count}.IndexJ = dof+IndexJ;

A{count}.A1 = A1;%-G(row1,col1);
A{count}.A2 = A2;% GT(row2,col2);%G(row2,col2)';
A{count}.A3 = A3;%I3(row3,col3);%
% end

A1 = H(row1,col1);
A2 = GT(row2,col2);%G(row2,col2)';
A3 = I3(row3,col3);%
% if norm(A1-z1)>1e-4 && norm(A2-z2)>1e-4 && norm(A3-z3)>1e-4
count = count+1;
A{count}.IndexI = IndexI;
A{count}.IndexJ = dof+IndexJ;

A{count}.A1=A1;%H(row1,col1);
A{count}.A2=A2;%GT(row2,col2);%G(row2,col2)';
A{count}.A3=A3;%I3(row3,col3);%
% end

A1 = G(row1,col1);
A2 = K(row2,col2);%H';
A3 = I3(row3,col3);%
% if norm(A1-z1)>1e-4 && norm(A2-z2)>1e-4 && norm(A3-z3)>1e-4
count = count+1;
A{count}.IndexI = IndexI;
A{count}.IndexJ = dof+IndexJ;

A{count}.A1=A1;%G(row1,col1);
A{count}.A2=A2;%K(row2,col2);%H';
A{count}.A3=A3;%I3(row3,col3);%
% end

% A13
A1 = -G(row1,col1);
A2 = I2(row2,col2);%
A3 = GT(row3,col3);%G(row3,col3)';
% if norm(A1-z1)>1e-4 && norm(A2-z2)>1e-4 && norm(A3-z3)>1e-4
count = count+1;
A{count}.IndexI = IndexI;
A{count}.IndexJ = 2*dof+IndexJ;

A{count}.A1=A1;%-G(row1,col1);
A{count}.A2=A2;%I2(row2,col2);%
A{count}.A3=A3;%GT(row3,col3);%G(row3,col3)';
% end

A1 = H(row1,col1);
A2 = I2(row2,col2);%
A3 = GT(row3,col3);%G(row3,col3)';
% if norm(A1-z1)>1e-4 && norm(A2-z2)>1e-4 && norm(A3-z3)>1e-4
count = count+1;
A{count}.IndexI = IndexI;
A{count}.IndexJ = 2*dof+IndexJ;

A{count}.A1=A1;%H(row1,col1);
A{count}.A2=A2;%I2(row2,col2);%
A{count}.A3=A3;%GT(row3,col3);%G(row3,col3)';
% end

A1 = G(row1,col1);
A2 = I2(row2,col2);%
A3 = K(row3,col3);%H';
% if norm(A1-z1)>1e-4 && norm(A2-z2)>1e-4 && norm(A3-z3)>1e-4
count = count+1;
A{count}.IndexI = IndexI;
A{count}.IndexJ = 2*dof+IndexJ;

A{count}.A1=A1;%G(row1,col1);
A{count}.A2=A2;%I2(row2,col2);%
A{count}.A3=A3;%K(row3,col3);%H';
% end

% A21
A1 = -GT(row1,col1);%G(row1,col1)';
A2 = G(row2,col2);
A3 = I3(row3,col3);%
% if norm(A1-z1)>1e-4 && norm(A2-z2)>1e-4 && norm(A3-z3)>1e-4
count = count+1;
A{count}.IndexI = dof+IndexI;
A{count}.IndexJ = IndexJ;

A{count}.A1=A1;%-GT(row1,col1);%G(row1,col1)';
A{count}.A2=A2;%G(row2,col2);
A{count}.A3=A3;%I3(row3,col3);%
% end

A1 = GT(row1,col1);%G(row1,col1)';
A2 = H(row2,col2);
A3 = I3(row3,col3);%
% if norm(A1-z1)>1e-4 && norm(A2-z2)>1e-4 && norm(A3-z3)>1e-4
count = count+1;
A{count}.IndexI = dof+IndexI;
A{count}.IndexJ = IndexJ;

A{count}.A1=A1;%GT(row1,col1);%G(row1,col1)';
A{count}.A2=A2;%H(row2,col2);
A{count}.A3=A3;%I3(row3,col3);%
% end

A1 = K(row1,col1);%H';
A2 = G(row2,col2);
A3 = I3(row3,col3);%
% if norm(A1-z1)>1e-4 && norm(A2-z2)>1e-4 && norm(A3-z3)>1e-4
count = count+1;
A{count}.IndexI = dof+IndexI;
A{count}.IndexJ = IndexJ;

A{count}.A1=A1;%K(row1,col1);%H';
A{count}.A2=A2;%G(row2,col2);
A{count}.A3=A3;%I3(row3,col3);%
% end
% A22
A1 = GG(row1,col1)-J(row1,col1)-L(row1,col1)+alpha*Q(row1,col1)-w2*I1(row1,col1);%speye(nx1,ny1);
A2 = I2(row2,col2);%
A3 = I3(row3,col3);%
% if norm(A1-z1)>1e-4 && norm(A2-z2)>1e-4 && norm(A3-z3)>1e-4
count = count+1;
A{count}.IndexI = dof+IndexI;
A{count}.IndexJ = dof+IndexJ;

A{count}.A1=A1;%GG(row1,col1)-J(row1,col1)-L(row1,col1)+alpha*Q(row1,col1)-w2*I1(row1,col1);%speye(nx1,ny1);
A{count}.A2=A2;%I2(row2,col2);%
A{count}.A3=A3;%I3(row3,col3);%
% end

A1 = I1(row1,col1);%
A2 = I2(row2,col2);%
A3 = GG(row3,col3)-J(row3,col3)-L(row3,col3)+alpha*Q(row3,col3);
% if norm(A1-z1)>1e-4 && norm(A2-z2)>1e-4 && norm(A3-z3)>1e-4
count = count+1;
A{count}.IndexI = dof+IndexI;
A{count}.IndexJ = dof+IndexJ;

A{count}.A1=A1;%I1(row1,col1);%
A{count}.A2=A2;%I2(row2,col2);%
A{count}.A3=A3;%GG(row3,col3)-J(row3,col3)-L(row3,col3)+alpha*Q(row3,col3);
% end

% A23
A1 = -I1(row1,col1);%
A2 = G(row2,col2);
A3 = GT(row3,col3);%G(row3,col3)';
% if norm(A1-z1)>1e-4 && norm(A2-z2)>1e-4 && norm(A3-z3)>1e-4
count = count+1;
A{count}.IndexI = dof+IndexI;
A{count}.IndexJ = 2*dof+IndexJ;

A{count}.A1=A1;%-I1(row1,col1);%
A{count}.A2=A2;%G(row2,col2);
A{count}.A3=A3;%GT(row3,col3);%G(row3,col3)';
% end

A1 = I1(row1,col1);%
A2 = H(row2,col2);
A3 = GT(row3,col3);%G(row3,col3)';
% if norm(A1-z1)>1e-4 && norm(A2-z2)>1e-4 && norm(A3-z3)>1e-4
count = count+1;
A{count}.IndexI = dof+IndexI;
A{count}.IndexJ = 2*dof+IndexJ;

A{count}.A1=A1;%I1(row1,col1);%
A{count}.A2=A2;%H(row2,col2);
A{count}.A3=A3;%GT(row3,col3);%G(row3,col3)';
% end

A1 = I1(row1,col1);%
A2 = G(row2,col2);
A3 = K(row3,col3);%H';
% if norm(A1-z1)>1e-4 && norm(A2-z2)>1e-4 && norm(A3-z3)>1e-4
count = count+1;
A{count}.IndexI = dof+IndexI;
A{count}.IndexJ = 2*dof+IndexJ;

A{count}.A1=A1;%I1(row1,col1);%
A{count}.A2=A2;%G(row2,col2);
A{count}.A3=A3;%K(row3,col3);%H';
% end

% A31
A1 = -GT(row1,col1);%G(row1,col1)';
A2 = I2(row2,col2);%
A3 = G(row3,col3);
% if norm(A1-z1)>1e-4 && norm(A2-z2)>1e-4 && norm(A3-z3)>1e-4
count = count+1;
A{count}.IndexI = 2*dof+IndexI;
A{count}.IndexJ = IndexJ;

A{count}.A1=A1;%-GT(row1,col1);%G(row1,col1)';
A{count}.A2=A2;%I2(row2,col2);%
A{count}.A3=A3;%G(row3,col3);
% end

A1 = GT(row1,col1);%G(row1,col1)';
A2 = I2(row2,col2);%
A3 = H(row3,col3);
% if norm(A1-z1)>1e-4 && norm(A2-z2)>1e-4 && norm(A3-z3)>1e-4
count = count+1;
A{count}.IndexI = 2*dof+IndexI;
A{count}.IndexJ = IndexJ;

A{count}.A1=A1;%GT(row1,col1);%G(row1,col1)';
A{count}.A2=A2;%I2(row2,col2);%
A{count}.A3=A3;%H(row3,col3);
% end

A1 = K(row1,col1);%H';
A2 = I2(row2,col2);%
A3 = G(row3,col3);
% if norm(A1-z1)>1e-4 && norm(A2-z2)>1e-4 && norm(A3-z3)>1e-4
count = count+1;
A{count}.IndexI = 2*dof+IndexI;
A{count}.IndexJ = IndexJ;

A{count}.A1=A1;%K(row1,col1);%H';
A{count}.A2=A2;%I2(row2,col2);%
A{count}.A3=A3;%G(row3,col3);
% end

% A32
A1 = -I1(row1,col1);%
A2 = GT(row2,col2);%G(row2,col2)';
A3 = G(row3,col3);
% if norm(A1-z1)>1e-4 && norm(A2-z2)>1e-4 && norm(A3-z3)>1e-4
count = count+1;
A{count}.IndexI = 2*dof+IndexI;
A{count}.IndexJ = dof+IndexJ;

A{count}.A1=A1;%-I1(row1,col1);%
A{count}.A2=A2;%GT(row2,col2);%G(row2,col2)';
A{count}.A3=A3;%G(row3,col3);
% end

A1 = I1(row1,col1);%
A2 = GT(row2,col2);%G(row2,col2)';
A3 = H(row3,col3);
% if norm(A1-z1)>1e-4 && norm(A2-z2)>1e-4 && norm(A3-z3)>1e-4
count = count+1;
A{count}.IndexI = 2*dof+IndexI;
A{count}.IndexJ = dof+IndexJ;

A{count}.A1=A1;%I1(row1,col1);%
A{count}.A2=A2;%GT(row2,col2);%G(row2,col2)';
A{count}.A3=A3;%H(row3,col3);
% end

A1 = I1(row1,col1);%
A2 = K(row2,col2);%H';
A3 = G(row3,col3);
% if norm(A1-z1)>1e-4 && norm(A2-z2)>1e-4 && norm(A3-z3)>1e-4
count = count+1;
A{count}.IndexI = 2*dof+IndexI;
A{count}.IndexJ = dof+IndexJ;

A{count}.A1=A1;%I1(row1,col1);%
A{count}.A2=A2;%K(row2,col2);%H';
A{count}.A3=A3;%G(row3,col3);
% end

% A33
A1 = GG(row1,col1)-J(row1,col1)-L(row1,col1)+alpha*Q(row1,col1)-w2*I1(row1,col1);%speye(nx1,ny1);
A2 = I2(row2,col2);%
A3 = I3(row3,col3);%
% if norm(A1-z1)>1e-4 && norm(A2-z2)>1e-4 && norm(A3-z3)>1e-4
count = count+1;
A{count}.IndexI = 2*dof+IndexI;
A{count}.IndexJ = 2*dof+IndexJ;

A{count}.A1=A1;%GG(row1,col1)-J(row1,col1)-L(row1,col1)+alpha*Q(row1,col1)-w2*I1(row1,col1);%speye(nx1,ny1);
A{count}.A2=A2;%I2(row2,col2);%
A{count}.A3=A3;%I3(row3,col3);%
% end

A1 = I1(row1,col1);%
A2 = GG(row2,col2)-J(row2,col2)-L(row2,col2)+alpha*Q(row2,col2);
A3 = I3(row3,col3);%
% if norm(A1-z1)>1e-4 && norm(A2-z2)>1e-4 && norm(A3-z3)>1e-4
count = count+1;
A{count}.IndexI = 2*dof+IndexI;
A{count}.IndexJ = 2*dof+IndexJ;

A{count}.A1=A1;%I1(row1,col1);%
A{count}.A2=A2;%GG(row2,col2)-J(row2,col2)-L(row2,col2)+alpha*Q(row2,col2);
A{count}.A3=A3;%I3(row3,col3);%
% end