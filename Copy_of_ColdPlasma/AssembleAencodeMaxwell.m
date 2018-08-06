function A_encode = AssembleAencodeMaxwell(IndexI,IndexJ,GG,J,H,K,L,Q,II) 
% Assemble the global matrix
% 

   % A11
    count = 1;
    A_encode{count}.IndexI = IndexI';
    A_encode{count}.IndexJ = IndexJ';
    
    A_encode{count}.A1=II;
    A_encode{count}.A2=GG-c1*J-c2*L+alpha*Q-pde.w2*II;
    A_encode{count}.A3=II;
%     full(sparse(IndexI',IndexJ',kron(II,kron(A_encode{count}.A2,II))))
    
    count = count+1;
    A_encode{count}.IndexI = IndexI';
    A_encode{count}.IndexJ = IndexJ';
    
    A_encode{count}.A1=II;
    A_encode{count}.A2=II;
    A_encode{count}.A3=GG-c1*J-c2*L+alpha*Q;
    
%     full(sparse(IndexI',IndexJ',kron(II,kron(II,A_encode{count}.A3))))
    
    % A12
    count = count+1;
    A_encode{count}.IndexI = IndexI';
    A_encode{count}.IndexJ = dof_1D_x^3+IndexJ';
    
    A_encode{count}.A1=-G;
    A_encode{count}.A2= G';
    A_encode{count}.A3=II;
    
    count = count+1;
    A_encode{count}.IndexI = IndexI';
    A_encode{count}.IndexJ = dof_1D_x^3+IndexJ';
    
    A_encode{count}.A1=H;
    A_encode{count}.A2=G';
    A_encode{count}.A3=c1*II;
    
    count = count+1;
    A_encode{count}.IndexI = IndexI';
    A_encode{count}.IndexJ = dof_1D_x^3+IndexJ';
    
    A_encode{count}.A1=G;
    A_encode{count}.A2=K;%H';
    A_encode{count}.A3=c2*II;
    
    % A13
    count = count+1;
    A_encode{count}.IndexI = IndexI';
    A_encode{count}.IndexJ = 2*dof_1D_x^3+IndexJ';
    
    A_encode{count}.A1=-G;
    A_encode{count}.A2=II;
    A_encode{count}.A3=G';
    
    count = count+1;
    A_encode{count}.IndexI = IndexI';
    A_encode{count}.IndexJ = 2*dof_1D_x^3+IndexJ';
    
    A_encode{count}.A1=H;
    A_encode{count}.A2=c1*II;
    A_encode{count}.A3=G';
    
    count = count+1;
    A_encode{count}.IndexI = IndexI';
    A_encode{count}.IndexJ = 2*dof_1D_x^3+IndexJ';
    
    A_encode{count}.A1=G;
    A_encode{count}.A2=c2*II;
    A_encode{count}.A3=K;%H';
    
    % A21
    count = count+1;
    A_encode{count}.IndexI = dof_1D_x^3+IndexI';
    A_encode{count}.IndexJ = IndexJ';
    
    A_encode{count}.A1=-G';
    A_encode{count}.A2=G;
    A_encode{count}.A3=II;
    
    count = count+1;
    A_encode{count}.IndexI = dof_1D_x^3+IndexI';
    A_encode{count}.IndexJ = IndexJ';
    
    A_encode{count}.A1=G';
    A_encode{count}.A2=H;
    A_encode{count}.A3=c1*II;
    
    count = count+1;
    A_encode{count}.IndexI = dof_1D_x^3+IndexI';
    A_encode{count}.IndexJ = IndexJ';
    
    A_encode{count}.A1=K;%H';
    A_encode{count}.A2=G;
    A_encode{count}.A3=c2*II;
    
    % A22
    count = count+1;
    A_encode{count}.IndexI = dof_1D_x^3+IndexI';
    A_encode{count}.IndexJ = dof_1D_x^3+IndexJ';
    
    A_encode{count}.A1=GG-c1*J-c2*L+alpha*Q-pde.w2*II;
    A_encode{count}.A2=II;
    A_encode{count}.A3=II;
    
    count = count+1;
    A_encode{count}.IndexI = dof_1D_x^3+IndexI';
    A_encode{count}.IndexJ = dof_1D_x^3+IndexJ';
    
    A_encode{count}.A1=II;
    A_encode{count}.A2=II;
    A_encode{count}.A3=GG-c1*J-c2*L+alpha*Q;
    
    % A23
    count = count+1;
    A_encode{count}.IndexI = dof_1D_x^3+IndexI';
    A_encode{count}.IndexJ = 2*dof_1D_x^3+IndexJ';
    
    A_encode{count}.A1=-II;
    A_encode{count}.A2=G;
    A_encode{count}.A3=G';
    
    count = count+1;
    A_encode{count}.IndexI = dof_1D_x^3+IndexI';
    A_encode{count}.IndexJ = 2*dof_1D_x^3+IndexJ';
    
    A_encode{count}.A1=c1*II;
    A_encode{count}.A2=H;
    A_encode{count}.A3=G';
    
    count = count+1;
    A_encode{count}.IndexI = dof_1D_x^3+IndexI';
    A_encode{count}.IndexJ = 2*dof_1D_x^3+IndexJ';
    
    A_encode{count}.A1=c2*II;
    A_encode{count}.A2=G;
    A_encode{count}.A3=K;%H';
    
    % A31
    count = count+1;
    A_encode{count}.IndexI = 2*dof_1D_x^3+IndexI';
    A_encode{count}.IndexJ = IndexJ';
    
    A_encode{count}.A1=-G';
    A_encode{count}.A2=II;
    A_encode{count}.A3=G;
    
    count = count+1;
    A_encode{count}.IndexI = 2*dof_1D_x^3+IndexI';
    A_encode{count}.IndexJ = IndexJ';
    
    A_encode{count}.A1=G';
    A_encode{count}.A2=c1*II;
    A_encode{count}.A3=H;
    
    count = count+1;
    A_encode{count}.IndexI = 2*dof_1D_x^3+IndexI';
    A_encode{count}.IndexJ = IndexJ';
    
    A_encode{count}.A1=K;%H';
    A_encode{count}.A2=c2*II;
    A_encode{count}.A3=G;
    
    % A32
    count = count+1;
    A_encode{count}.IndexI = 2*dof_1D_x^3+IndexI';
    A_encode{count}.IndexJ = dof_1D_x^3+IndexJ';
    
    A_encode{count}.A1=-II;
    A_encode{count}.A2=G';
    A_encode{count}.A3=G;
    
    count = count+1;
    A_encode{count}.IndexI = 2*dof_1D_x^3+IndexI';
    A_encode{count}.IndexJ = dof_1D_x^3+IndexJ';
    
    A_encode{count}.A1=c1*II;
    A_encode{count}.A2=G';
    A_encode{count}.A3=H;
    
    count = count+1;
    A_encode{count}.IndexI = 2*dof_1D_x^3+IndexI';
    A_encode{count}.IndexJ = dof_1D_x^3+IndexJ';
    
    A_encode{count}.A1=c2*II;
    A_encode{count}.A2=K;%H';
    A_encode{count}.A3=G;
    
    % A33
    count = count+1;
    A_encode{count}.IndexI = 2*dof_1D_x^3+IndexI';
    A_encode{count}.IndexJ = 2*dof_1D_x^3+IndexJ';
    
    A_encode{count}.A1=GG-c1*J-c2*L+alpha*Q-pde.w2*II;
    A_encode{count}.A2=II;
    A_encode{count}.A3=II;
    
    count = count+1;
    A_encode{count}.IndexI = 2*dof_1D_x^3+IndexI';
    A_encode{count}.IndexJ = 2*dof_1D_x^3+IndexJ';
    
    A_encode{count}.A1=II;
    A_encode{count}.A2=GG-c1*J-c2*L+alpha*Q;
    A_encode{count}.A3=II;