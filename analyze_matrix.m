function analyze_matrix(A)

fprintf('\n\nbegin matrix analysis\n');

% these tests disabled - they are slow, and none of the problem
% matrices tested are symmetric or positive definite

% 1: test for symmetry
%is_sym = isequal(A,A');


%if is_sym
%    fprintf("--symmetric--\n");
%else
%    fprintf("--nonsymmetric--\n");
%end



% 2: test for pos def
%test = zeros(size(A));
%if is_sym
%    test = A;
%else
%    test = A + A';
%end

%pos_def = 0;
%eigenvals = eig(test);
%non_positive_eigenvals = (eigenvals <= 0);
%if sum(non_positive_eigenvals) > 0
%    pos_def = 1;
%end

%if pos_def    
%    fprintf("--positive definite--\n");
%else
%    fprintf("--not positive definite--\n");
%end

 
% 3: spectrum
[eigenvectors, eigenvalue_mat] = eig(A);
eigenvalues = diag(eigenvalue_mat);
figure
plot(eigenvalues,'o','MarkerSize',12);
saveas(gcf,'spectrum.png');


% 4: condition number of eigenvect mat
fprintf('eigen condition number: %d\n', rcond(eigenvectors));

% 5: condition number of matrix
fprintf('matrix condition number: %d\n', rcond(A));

% 6: sparsity %
fprintf('matrix density: %f\n', nnz(A)/numel(A));

% 7: sparsity pattern
spy(A);

fprintf('end matrix analysis\n\n');

end