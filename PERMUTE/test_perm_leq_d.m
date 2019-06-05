
disp('simple test of perm_leq_d ');
nerrors = 0;

idebug = 0;

idim = 3;
LevMax = [2 3 4];
SumMax = 4;
disp(sprintf('sum of indices of dimension %d <=  %d', ...
     idim,  SumMax ));

for icase=1:3,
for increasing_sum_order=0:1,
   if (icase == 1),
       LevMax_in = min(LevMax);
   end;
   if (icase == 2),
       LevMax_in = max(LevMax);
   end;
   if (icase == 3),
       LevMax_in = LevMax;
   end;

   disp(sprintf('increasing_sum_order = %d', increasing_sum_order ));


   icount = perm_leq_d_count(idim,LevMax_in,SumMax);
   result = perm_leq_d(idim,LevMax_in,SumMax, increasing_sum_order)

   if (idebug >= 1),
   
    disp(sprintf('perm_leq_d_count(idim=%d,SumMax=%d) returns icount=%d', ...
                                   idim, SumMax, icount ));
    disp(sprintf('size(result) = (%d,%d) ', ...
                  size(result,1), size(result,2) ));
   end;

   isok = (size(result,1) == icount);
   if (~isok),
           disp(sprintf('size(result,1)=%g, icount=%g', ...
                         size(result,1),  icount ));
           nerrors = nerrors + 1;
   end;

   is_scalar = (size(LevMax_in,1)==1) && (size(LevMax_in,2)==1);
   isok_dim = zeros(idim,1);
   for i=1:idim,
      lmax = LevMax(i);
      if (is_scalar),
         lmax = LevMax_in;
      end;
      isok_dim(i) = all( (0 <= result(:,i)) & ...
                      (result(:,i) <= lmax) );
   end;

   isok = all(isok_dim(1:idim)) && ...
          all(sum( result,2) <= SumMax);
   if (~isok),
       idx = find( sum( result,2) > SumMax);
       for i=1:numel(idx),
               j = idx(i);
               disp(sprintf('j=%d,result(j,:) ',j));
               result(j, 1:idim)
       end;
       nerrors = nerrors + 1;
   end;

end;
end;
if (nerrors == 0),
        disp('ALL OK');
else
   disp(sprintf('There are %d errors ',nerrors));
end;
