disp('simple test of perm_eq ');
idim = 3;
n = 4;
disp(sprintf('sum of indices of dimension %d equals %d', ...
     idim,  n ));

result = perm_eq(idim,n)
icount = perm_eq_count(idim,n);

disp(sprintf('perm_eq_count(idim=%d,n=%d) returns icount=%d', ...
                      idim, n, icount ));
disp(sprintf('size(result) = (%d,%d) ', ...
       size(result,1), size(result,2) ));
