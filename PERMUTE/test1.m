disp('simple test of perm_eq ');
nerrors = 0;
idim = 3;
n = 4;
disp(sprintf('sum of indices of dimension %d equals %d', ...
     idim,  n ));

for last_index_decreasing=0:1,
  disp(sprintf('last_index_decreasing=%d', ...
                last_index_decreasing ));

  icount = perm_eq_count(idim,n);
  result = perm_eq(idim,n, last_index_decreasing)
  
  disp(sprintf('perm_eq_count(idim=%d,n=%d) returns icount=%d', ...
                      idim, n, icount ));
  disp(sprintf('size(result) = (%d,%d) ', ...
       size(result,1), size(result,2) ));

  isok = (size(result,1) == icount);
  if (~isok),
          disp(sprintf('size(result,1)=%d, icount=%d',size(result,1), icount ));
          nerrors = nerrors + 1;
  end;

  isok = all(sum( result, 2) == n);
  if (~isok),
        idx = find(sum(result,2) ~= n)
        for i=1:numel(idx),
            j = idx(i);
            disp(sprintf('j=%d, result(j,:) = ', j));
            result( j, 1:idim )
        end;
        nerrors = nerrors + 1;
  end;
end;

if (nerrors == 0),
        disp('ALL OK');
end;
