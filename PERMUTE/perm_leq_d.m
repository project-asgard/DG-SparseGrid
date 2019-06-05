function [result] = perm_leq_d( idim, LevMax, SumMax, order_by_n_in )
% [result] = perm_leq( idim, n  [, order_by_n_in] )
%
% return tuples where sum of indices is less than or
% equal to n
%
icount = perm_leq_count( idim, n );
result = zeros( icount, idim );

order_by_n = 0;
if (nargin >= 3),
    order_by_n = order_by_n_in;
end;


if (order_by_n),
    ip = 1;
    for i=0:SumMax,
        SumEq = i;
        isize = perm_eq_d_count(idim, LevMax,SumEq);
        
        result( ip:(ip+isize-1),1:idim) = perm_eq_d(idim,LevMax,i);
        
        ip = ip + isize;
    end;
else
    if (idim == 1),
        n = LevMax(1);
        result = reshape(0:n, (n+1),1);
        return;
    end;
    
    n = LevMax(idim); 
    ip = 1;
    for i=0:min(n,SumMax),
        isize = perm_leq_d_count(idim-1, n-i);
        result( ip:(ip+isize-1),1:(idim-1)) = ...
            perm_leq_d(idim-1,LevMax(1:(idim-1))-i,SumMax-i,order_by_n);
        result( ip:(ip+isize-1),idim) = i;
        
        ip = ip + isize;
    end;
end;

return
end

