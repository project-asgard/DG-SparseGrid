% The pattern of Sparsity

figure(3)
subplot(1,2,1)
spy(A_2D)
title('A Sparse Full-Grid Matrix')
nz = nnz(A_2D);
xlabel(sprintf('nonzeros = %d (%.3f%%)',nz,nz*100/numel(A_2D)));
subplot(1,2,2)
spy(A_s)
title('A Sparse Sparse-Grid Matrix')
nz = nnz(A_s);
xlabel(sprintf('nonzeros = %d (%.3f%%)',nz,nz*100/numel(A_s)));