function result = index_max( ndim, Nmat, Levmax)
% result = index_max( ndim, Nmat, Levmax)
%
% Nmat is cell array of size ndim
% Nmat{i}  is a list of integers
%
% if ndim == 4,
%
% result(i,1:4) is such that  
% Nmat{1}(i1) + Nmat{2}(i2) + ... Nmat{4}(i4)  <= Levmax
% i1 = result(i,1); i2 = result(i,2); i3 = result(i,3); i4 = result(i,4);
%
if (ndim == 1),
        result = find( Nmat{1} <= Levmax );
        result = reshape( result, numel(result),1 );
        return;
else
        m = index_max_count( ndim, Nmat, Levmax );
        result = zeros( m, ndim);

        v = Nmat{ndim};
        ip = 1;

        is_valid = find( v <= Levmax );
        isize = index_max_count( ndim-1,Nmat,Levmax);
        result1 = index_max( ndim-1,Nmat,Levmax);
        for i=1:numel(is_valid),
            
          result(ip:(ip+isize-1), 1:(ndim-1)) = ...
                        result1(1:isize,1:(ndim-1));

          vi = is_valid(i);
          result(ip:(ip+isize-1),ndim) = vi;

          ip = ip + isize;
         end;
end

end
        


