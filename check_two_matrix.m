close all
Dof = Deg^2*size(HASHInv,2);
Dim = 2;
AA_new = sparse(Dof,Dof);
for i = 1:size(A_new,2)

    A = A_new{i}.A{1};
    B = A_new{i}.A{2};
    II =  A_new{i}.IndexI;
    JJ =  A_new{i}.IndexJ;
    
    AA_new(II,JJ) = AA_new(II,JJ) +kron(A,B);
end

AA_old = sparse(Dof,Dof);
for i = 1:size(A_old,2)
    A = A_old{i}.A1;
    B = A_old{i}.A2;
    II =  A_old{i}.IndexI;
    JJ =  A_old{i}.IndexJ;
    
    AA_old(II,JJ) = AA_old(II,JJ) +kron(A,B);
end
    
figure;
subplot(1,3,1); spy(AA_new);title('New')
subplot(1,3,2); spy(AA_old);title('Old')
B = AA_old-AA_new;
B = abs(B)>1e-10;
subplot(1,3,3); spy(B);

b = ones(Dof,1);
figure
plot(AA_old*b-AA_new*b)