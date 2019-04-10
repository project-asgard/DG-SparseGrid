function m = index_max_count(ndim, Nmat, Levmax)
% m = index_max_count(ndim, Nmat, Levmax)
%
% just count number of entries that would be returned
% by index_max()
%
if (ndim == 1),
        m = numel( find( Nmat{1} <= Levmax ) );
else
        m = 0;
        is_valid = find( Nmat{ndim} <= Levmax );
        if (numel(is_valid) >= 1),
          mtemp = index_max_count(ndim-1,Nmat,Levmax);
          m = numel( is_valid) * mtemp;
        end;
end

end

