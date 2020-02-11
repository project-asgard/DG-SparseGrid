function [x, iter, relres, relres_hist] = iterative_solve(pde, opts, A, b, tol, max_iter, restart, graph_residual, precond)
%iterative solve driver

% set options
if ~exist('restart','var')
    restart = [];
else
   if strcmpi(opts.solve_choice, 'BICGSTAB')
       disp('non-restarting method')
   end
end

if ~exist('tol','var') || isempty(tol)
    tol = 1e-6; % matlab default
end

if ~exist('max_iter', 'var') || isempty(max_iter)
    max_iter = size(A, 1);
end

% set up preconditioner
if strcmpi(opts.preconditioner, 'JACOBI') %point jacobi
    M = diag(diag(A));
elseif strcmpi(opts.preconditioner, 'BLK_JACOBI') % block jacobi
   blk_n = pde.deg^size(pde.dimensions,1);
   nblks = size(A,1) / blk_n; 
   kron_map = kron(eye(nblks), ones(blk_n, blk_n));
   blk_list = A(kron_map ~= 0);
   blocks = reshape(blk_list,blk_n,blk_n,size(A,1)/blk_n);
   M = num2cell(blocks, [1,2]);
   M = reshape(M, 1, size(M, 3));
   M = blkdiag(M{:});
elseif strcmpi(opts.preconditioner, 'CUSTOM')
    if ~exist('precond', 'var')
        error('selected custom preconditioning, did not supply one');
    end
    M = precond; %custom preconditioner
else
    M = [];
end

% do solve
tic;
if strcmpi(opts.solve_choice, 'GMRES')
    [x,flag, relres, iter, relres_hist] = gmres(A, b, restart, tol, max_iter, M);
    if flag == 0; fprintf('GMRES converged w iters\n'); disp(iter);
    else; fprintf('GMRES failed with flag %d\n', flag); end
elseif strcmpi(opts.solve_choice, 'BICGSTAB')
    [x,flag, relres, iter, relres_hist] = bicgstab(A, b, tol, max_iter, M);
    if flag == 0; fprintf('BIGCGSTAB converged w iters\n'); disp(iter);
    else; fprintf('BICGSTAB failed with flag %d\n', flag); end
else
    error('solver choice not recongized!');
end
disp(['solve duration: ', num2str(toc), 's']);


% graph performance
if graph_residual
    figure(2000);
    semilogy(1:size(relres_hist,1),relres_hist/norm(b),'-o');
    xlabel('Iteration number');
    ylabel('Relative residual');
    
end
    
end

