function [result] = perm_leq_d( idim, LevMax_in, SumMax, ...
                                increasing_sum_order_in )
% [result] = perm_leq_d( idim, LevMax, SumMax, ...
%                        [, increasing_sum_order_in] )
%
% return tuples where sum of indices is less than or
% equal to SumMax and  i-th index in [0,LevMax(i)]
%

LevMax = LevMax_in;
is_scalar = (size(LevMax_in,1)==1) && (size(LevMax_in,2)==1);
if (is_scalar),
    LevMax  = LevMax_in * ones(idim,1);
end;


icount = perm_leq_d_count( idim, LevMax, SumMax );
result = zeros( icount, idim );

increasing_sum_order = 0;
if (nargin >= 4),
    increasing_sum_order = increasing_sum_order_in;
end;


if (increasing_sum_order),
    i1 = 0;
    inc = 1;
    i2 = SumMax;
else
    i1 = SumMax;
    inc = -1;
    i2 = 0;
end;

    ip = 1;
    for i=i1:inc:i2,
        SumEq = i;
        isize = perm_eq_d_count(idim, LevMax,SumEq);
        if (isize >= 1),
         result( ip:(ip+isize-1),1:idim) = perm_eq_d(idim,LevMax, SumEq);
         ip = ip + isize;
        end;
    end;
return
end

