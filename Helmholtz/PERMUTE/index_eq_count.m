function m = index_eq_count(ndim, Nmat, Levsum)
% m = index_eq_count(ndim, Nmat, Levsum)
%
% just count number of entries that would be returned
% by index_eq()
%
if (ndim == 1),
        m = numel( find( Nmat{1} == Levsum ) );
else
        m = 0;
        v = Nmat{ndim};
        for i=1:numel(v),
                mtemp = index_eq_count(ndim-1,Nmat,Levsum - v(i));
                m = m + mtemp;
        end;
end

end

