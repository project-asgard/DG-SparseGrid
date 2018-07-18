
disp('simple test of perm_leq ');
idim = 3;
n = 4;
disp(sprintf('sum of indices of dimension %d <=  %d', ...
     idim,  n ));

result = perm_leq(idim,n)
icount = perm_leq_count(idim,n);

disp(sprintf('perm_leq_count(idim=%d,n=%d) returns icount=%d', ...
                      idim, n, icount ));
disp(sprintf('size(result) = (%d,%d) ', ...
       size(result,1), size(result,2) ));
