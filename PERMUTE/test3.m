
disp('simple test of perm_max ');
idim = 3;
n = 4;
disp(sprintf('max of indices of dimension %d is  %d', ...
     idim,  n ));

result = perm_max(idim,n)
icount = perm_max_count(idim,n);

disp(sprintf('perm_max_count(idim=%d,n=%d) returns icount=%d', ...
                      idim, n, icount ));
disp(sprintf('size(result) = (%d,%d) ', ...
       size(result,1), size(result,2) ));
