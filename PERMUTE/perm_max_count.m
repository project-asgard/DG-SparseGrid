function icount = perm_max_count(idim,nvec)
%
% icount = perm_max_count(idim,n)
%
% number of tuples for dimension dim
% and max value of index <= nvec(1:idim)
%
is_scalar = (size(nvec,1)==1) && (size(nvec,2)==1);
if (is_scalar),
  icount = (nvec+1)^idim;
else
  icount = prod( nvec(1:idim) + 1 );
end;

return;
end

