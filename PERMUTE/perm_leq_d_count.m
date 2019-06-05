function [icount] = perm_leq_d_count( idim, LevMax, SumMax )
%
% [icount] = perm_leq_d_count( idim, LevMax, SumMax )
%
% compute number of tuples where
% dimension is idim and sum of indices is <= SumMax
% and i-th index in [0,LevMax(i)]
%
% for example idim = 2, LevMax = [2,2], SumMax = 2
%
% sum to 0:  (0,0)
% sum to 1:  (0,1), (1,0)
% sum to 2:  (0,2), (1,1), (2,0)
%
% thus icount is 1 + 2 + 3 = 6
%

icount = 0;
for lSumMax=0:SumMax,
  ires = perm_eq_d_count(idim, LevMax, lSumMax);
  icount = icount + ires;
end;

return;
end
