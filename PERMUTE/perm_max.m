function [result] = perm_max(idim,LevMax_in, last_index_decreasing_in)
%
% return tuples where max of indices  in i-th dimension <= LevMax(i)
%
% for example idim = 1, LevMax(1) = 2
%
% tuples in 1D are [0; 
%                   1; 
%                   2]
%
%
% for LevMax = [2,2], tuples in 2D are 
% 
% [0,0; 0,1; 0,2; 
%  1,0; 1,1; 1,2;
%  2,0; 2,1; 2,2]
%
LevMax = LevMax_in;
is_scalar = (size(LevMax,1) == 1) && (size(LevMax,2) == 1);
if (is_scalar)
        LevMax = ones(idim,1)*max(LevMax_in);
end;

last_index_decreasing = 0;
if (nargin >= 3),
        last_index_decreasing = last_index_decreasing_in;
end;

if (idim == 1),

 n = LevMax(idim);
 if (last_index_decreasing),
  result = reshape( n:-1:0, (n+1),1);
 else
  result = reshape(0:n, (n+1),1);
 end;

else;
 % --------------
 % here idim >= 2
 % --------------
  n = LevMax(idim);

  ivec = perm_max( idim-1,LevMax, last_index_decreasing );
  m = size(ivec,1);
  result = zeros( m*(n+1), idim );

  for i=0:n,
    i1 = i*m + 1;
    i2 = i1 + m - 1;
    if (last_index_decreasing),
      result( i1:i2, 1:(idim-1)) = ivec(1:m,1:(idim-1));
      result( i1:i2, idim) = n-i;
    else
      result( i1:i2, 1:(idim-1)) = ivec(1:m,1:(idim-1));
      result( i1:i2, idim) = i;
    end;
   end;
end;

return;

end

