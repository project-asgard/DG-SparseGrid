% approximation for full-grid of Poisson Eq

dof_full = dof^2;

A_2D = sparse(dof_full);
b_2D = sparse(dof_full,1);
sol_2D = sparse(dof_full,1);


% tic
A_2D = kron(M_mass,Stiff_1D)+kron(Stiff_1D,M_mass);
% toc
% full(A_2D)

b_2D = 2*pi^2*kron(b,b);

if dof_full<1e6
tic
sol_2D = A_2D\b_2D;
toc
else
    disp('Matirx is TOO BIG!')
    
end


cond_full = condest(A_2D);

Iu_2D = kron(M_mass,M_mass)\kron(uu,uu);

Lmax_2D=max(abs(Iu_2D-sol_2D));


