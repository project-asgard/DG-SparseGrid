function [result] = perm_max(idim,nvec_in, last_index_decreasing_in)
%
% return tuples where max of indices  in i-th dimension <= nvec(i)
%
% for example idim = 1, nvec(1) = 2
%
% tuples in 1D are [0; 
%                   1; 
%                   2]
%
%
% for nvec = [2,2], tuples in 2D are 
% 
% [0,0; 0,1; 0,2; 
%  1,0; 1,1; 1,2;
%  2,0; 2,1; 2,2]
%
nvec = nvec_in;
is_scalar = (size(nvec,1) == 1) && (size(nvec,2) == 1);
if (is_scalar)
        nvec = ones(idim,1)*max(nvec_in);
end;

last_index_decreasing = 0;
if (nargin >= 3),
        last_index_decreasing = last_index_decreasing_in;
end;

if (idim == 1),

 n = nvec(idim);
 if (last_index_decreasing),
  result = reshape( n:-1:0, (n+1),1);
 else
  result = reshape(0:n, (n+1),1);
 end;

else;
 % --------------
 % here idim >= 2
 % --------------
  n = nvec(idim);

  ivec = perm_max( idim-1,nvec, last_index_decreasing );
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

