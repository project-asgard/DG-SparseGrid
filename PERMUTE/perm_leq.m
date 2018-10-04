function [result] = perm_leq( idim, n )
%
% return tuples where sum of indices is less than or
% equal to n
%
icount = perm_leq_count( idim, n );
result = zeros( icount, idim );

ip = 1;
for i=0:n,
  isize = perm_eq_count(idim, i);

  result( ip:(ip+isize-1),:) = perm_eq(idim,i);

  ip = ip + isize;
end;

return
end
  
