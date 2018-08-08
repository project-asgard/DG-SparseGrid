function [result] = perm_eq(idim,n)
%
% return tuples where sum of indices equal to n
%
% for example idim = 2, n = 2
%
% tuples are (0,2), (1,1), (2,0)
% result is 3 by 2 matrix
% result = [0,2; ...
%           1,1; ...
%           2,0]
%

if (idim == 1),
  result = [n];
  return;
end;

if (idim == 2),
  ivec =  reshape(0:n,(n+1),1);
  result = zeros( n+1,2);
  result(:,1) = ivec;
  result(:,2) = n - ivec;
  return;
end;

icount = perm_eq_count(idim,n);
result = zeros( icount, idim);
ip = 1;
for i=0:n,
  isize = perm_eq_count( idim-1, i);
  result( ip:(ip+isize-1), 1:(idim-1)) = perm_eq( idim-1,i);
  result( ip:(ip+isize-1), idim ) = n-i;

  ip = ip + isize;
end;

return
end
