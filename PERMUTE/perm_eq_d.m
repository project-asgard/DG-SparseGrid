function [result] = perm_eq_d(idim,LevMax_in, SumEq, ...
                              last_index_decreasing_in)
% [result] = perm_eq_d(idim,LevMax, SumEq, last_index_decreasing_in)
%
% return tuples where sum of indices equal to SumEq
% and i-th index <= LevMax(i)
%
% for example idim = 2, LevMax = [3,3] SumEq = 2
%
% tuples are (0,2), (1,1), (2,0)
% result is 3 by 2 matrix
% result = [0,2; ...
%           1,1; ...
%           2,0]
%

LevMax = LevMax_in;
is_scalar = (size(LevMax_in,1) == 1) && (size(LevMax_in,2) == 1);
if (is_scalar),
    LevMax = LevMax_in * ones(idim,1);
end;

last_index_decreasing = 0;
if (nargin >= 4)
        last_index_decreasing = last_index_decreasing_in;
end;

if (idim == 1),
  if (LevMax(1) >= SumEq),
     result = [SumEq];
  else
     % ---------------------
     % impossible to satisfy
     % ---------------------
     result = [];
  end;
  return;
end;


icount = perm_eq_d_count(idim,LevMax,SumEq);
result = zeros( icount, idim);
ip = 1;

if (last_index_decreasing),
        i1 = min(LevMax(idim),SumEq);
        inc = -1;
        i2 = 0;
else
        i1 = 0;
        inc = 1;
        i2 = min(LevMax(idim),SumEq);
end;

for ilast=i1:inc:i2

  lSumEq = SumEq - ilast;

  isize = perm_eq_d_count( idim-1, LevMax(1:(idim-1)), lSumEq);

  tmp_result = perm_eq_d( idim-1,LevMax(1:(idim-1)), ...
                          lSumEq, last_index_decreasing );

  isok = (size(tmp_result,1) == isize);
  if (~isok),
     disp(sprintf('perm_eq_d: isize=%d,size(tmp_result,1)=%d', ...
                              isize,   size(tmp_result,1) ));
     disp(sprintf('idim-1=%d,lSumEq=%d',idim-1,lSumEq));
     disp('LevMax');
     LevMax(1:(idim-1))
     disp('tmp_result');
     tmp_result
  end;
  
  result( ip:(ip+isize-1), 1:(idim-1)) = tmp_result;

  result( ip:(ip+isize-1), idim ) = ilast;

  ip = ip + isize;
end;

return
end
