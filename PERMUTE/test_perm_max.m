
disp('simple test of perm_max ');
nerrors = 0;

idim = 3;
nvec = [2,3,4];
disp(sprintf('max of indices of dimension %d is ',idim));
nvec

for last_index_decreasing=0:1,
   disp(sprintf('last_index_decreasing=%d', last_index_decreasing));

   icount = perm_max_count(idim,nvec);
   result = perm_max(idim,nvec, last_index_decreasing)
  
   
   disp(sprintf('perm_max_count(idim=%d,nvec=[%g,%g,%g]) returns icount=%d', ...
     idim, nvec(1),nvec(2),nvec(3), icount ));
   disp(sprintf('size(result) = (%d,%d) ', ...
          size(result,1), size(result,2) ));
  
   isok = (size(result,1) == icount);
   if (~isok),
           disp(sprintf('size(result,1)=%g, icount=%d', ...
                         size(result,1),  icount ));

           nerrors = nerrors + 1;
   end;

  for i=1:idim
   n = nvec(i);
   isok = all(  result(:,i)  <= n );
   if (~isok),
           disp(sprintf('nvec(%d)=%g',i,nvec(i)));
           idx = find( result(:,i) > n );
           for i=1:numel(idx),
                   j = idx(i);
                   disp(sprintf('j=%d,result(j,:) ',j));
                   result(j,1:idim)
           end;
           nerrors = nerrors + 1;
   end;
  end;
end;

if (nerrors == 0),
        disp('ALL OK');
end;
