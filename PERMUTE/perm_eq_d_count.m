
function [icount] = perm_eq_d_count( idim, LevMax_in, SumEq )
%
% [icount] = perm_eq_d_count( idim, LevMax, SumEq )
%
% number of tuples for dimension = idim
% and sum of indices equal to SumEq
% subject to  i-th dimension index <= LevMax(i)
%
% for example dim = 2, LevMax = [3,3] SumEq = 3
% (0, 3), (1,2), (2, 1), (3,0)
% so icount is 4
%
LevMax = LevMax_in;
is_scalar = (size(LevMax_in,1) == 1) && (size(LevMax_in,2) == 1);
if (is_scalar),
        LevMax = LevMax_in * ones(idim,1);
end;



if (idim == 1)
  if (LevMax >= SumEq),
          icount = 1;
  else
          % ---------------------
          % impossible to satisfy
          % ---------------------
          icount = 0;
  end;

  return;
end;


icount = 0;
for ilast=0:min(LevMax(idim),SumEq),
  lSumEq = SumEq - ilast;
  ires = perm_eq_d_count( idim-1, LevMax(1:(idim-1)), lSumEq);
  icount = icount + ires;
end;

return
end
