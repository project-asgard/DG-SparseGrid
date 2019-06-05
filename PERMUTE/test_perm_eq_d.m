disp('simple test of perm_eq_d ');
nerrors = 0;
idim = 3;
SumEq = 5;
n = SumEq;
LevMax = [2,3,4];

disp('LevMax');
LevMax
idebug = 0;

disp(sprintf('sum of indices of dimension %d equals %d', ...
     idim,  SumEq ));

for last_index_decreasing=0:1,
  disp(sprintf('last_index_decreasing=%d', ...
                last_index_decreasing ));

  icount = perm_eq_d_count(idim,LevMax,SumEq);
  result = perm_eq_d(idim,LevMax, SumEq, last_index_decreasing)
  
  if (idebug >= 1),
   disp(sprintf('perm_eq_d_count(idim=%d,SumEq=%d) returns icount=%d', ...
                      idim, SumEq, icount ));
   disp(sprintf('size(result) = (%d,%d) ', ...
       size(result,1), size(result,2) ));
  end;

  isok = (size(result,1) == icount);
  if (~isok),
    disp(sprintf('size(result,1)=%d, icount=%d',...
                  size(result,1),    icount ));
    nerrors = nerrors + 1;
  end;

  isok_dim = zeros(idim,1);
  for i=1:idim,
     isok_dim(i) = all( (0 <= result(:,i)) & ...
                        (result(:,i) <= LevMax(i)));
  end;

  isok = all(isok_dim(1:idim)) && ...
         all(sum( result, 2) == SumEq);
         
  if (~isok),
        idx = find(sum(result,2) ~= SumEq)
        for i=1:numel(idx),
            j = idx(i);
            disp(sprintf('j=%d, result(j,:) = ', j));
            result( j, 1:idim )
        end;
        nerrors = nerrors + 1;

        disp(sprintf('idim=%d,SumEq=%d,last_index_decreasing=%d', ...
                      idim,   SumEq,   last_index_decreasing));
        disp('LevMax');
        LevMax(1:(idim))
  end;

end;

if (nerrors == 0),
        disp('ALL OK');
end;
