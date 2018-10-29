function [result] = perm_max(idim,n)
%
% return tuples where max of indices  <= n
%
% for example idim = 1, n = 2
%
% tuples in 1D are [0; 
%                   1; 
%                   2]
%
%
% tuples in 2D are 
% 
% [0,0; 0,1; 0,2; 
%  1,0; 1,1; 1,2;
%  2,0; 2,1; 2,2]
%

use_1st_index_fastest = 1;

if (idim == 1),
  result = reshape(0:n, (n+1),1);
else;
 % --------------
 % here idim >= 2
 % --------------

  ivec = perm_max( idim-1,n );
  m = size(ivec,1);
  result = zeros( m*(n+1), idim );

  for i=0:n,
    i1 = i*m + 1;
    i2 = i1 + m - 1;
    if (use_1st_index_fastest),
      result( i1:i2, 1:(idim-1)) = ivec;
      result( i1:i2, idim) = i;
    else
      result( i1:i2, 1) = i;
      result( i1:i2, 2:idim) = ivec;
    end;
   end;
end;

return;

end

