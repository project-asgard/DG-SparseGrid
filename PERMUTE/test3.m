
disp('simple test of perm_max ');
nerrors = 0;

idim = 3;
n = 4;
disp(sprintf('max of indices of dimension %d is  %d', ...
     idim,  n ));

for last_index_decreasing=0:1,
   icount = perm_max_count(idim,n);
   result = perm_max(idim,n, last_index_decreasing)
  
   
   disp(sprintf('perm_max_count(idim=%d,n=%d) returns icount=%d', ...
                         idim, n, icount ));
   disp(sprintf('size(result) = (%d,%d) ', ...
          size(result,1), size(result,2) ));
  
   isok = (size(result,1) == icount);
   if (~isok),
           disp(sprintf('size(result,1)=%g, icount=%d', ...
                         size(result,1),  icount ));

           nerrors = nerrors + 1;
   end;

   isok = all( max( result, [], 2) <= n );
   if (~isok),
           idx = find( max(result,[],2) > n );
           for i=1:numel(idx),
                   j = idx(i);
                   disp(sprintf('j=%d,result(j,:) ',j));
                   result(j,1:idim)
           end;
           nerrors = nerrors + 1;
   end;
end;

if (nerrors == 0),
        disp('ALL OK');
end;
