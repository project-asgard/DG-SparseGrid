
disp('simple test of perm_max ');
nerrors = 0;

idim = 3;
LevMax_org = [2,3,4];
disp(sprintf('max of indices of dimension %d is ',idim));
% LevMax

idebug = 0;
for icase=1:3,
for last_index_decreasing=0:1,
   if (icase == 1),
        LevMax = min(LevMax_org);
   end;
   if (icase == 2),
        LevMax = max(LevMax_org);
   end;
   if (icase == 3),
        LevMax = LevMax_org;
   end;

   disp(sprintf('last_index_decreasing=%d', last_index_decreasing));

   icount = perm_max_count(idim,LevMax);
   result = perm_max(idim,LevMax, last_index_decreasing)
  
   if (idebug >= 1), 
     disp(sprintf('perm_max_count(idim=%d) returns icount=%d', ...
                                  idim,  icount ));
     disp(sprintf('LevMax'));
     LevMax

     disp(sprintf('size(result) = (%d,%d) ', ...
                   size(result,1), size(result,2) ));
   end;
  
   isok = (size(result,1) == icount);
   if (~isok),
           disp(sprintf('size(result,1)=%g, icount=%d', ...
                         size(result,1),  icount ));

           nerrors = nerrors + 1;
   end;

  for i=1:idim
   is_scalar = (size(LevMax,1)==1) && (size(LevMax,2)==1);
   if (is_scalar),
       lmax = LevMax;
   else
       lmax = LevMax(i);
   end;

   isok = all(  (0 <= result(:,i)) & (result(:,i)  <= lmax ));
   if (~isok),
           disp(sprintf('LevMax(%d)=%g',i,LevMax(i)));
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
end;

if (nerrors == 0),
        disp('ALL OK');
else
        disp(sprintf('There are %d errors',nerrors));
end;
