function result = index_eq( ndim, Nmat, Levsum)
% result = index_eq( ndim, Nmat, Levsum)
%
% Nmat is cell array of size ndim
% Nmat{i}  is a list of integers
%
% if ndim == 4,
%
% result(i,1:4) is such that  
% Nmat{1}(i1) + Nmat{2}(i2) + ... Nmat{4}(i4)   = Levsum
% i1 = result(i,1); i2 = result(i,2); i3 = result(i,3); i4 = result(i,4);
%
if (ndim == 1),
        result = find( Nmat{1} == Levsum );
        result = reshape( result, numel(result),1 );
        return;
else
        m = index_eq_count( ndim, Nmat, Levsum );
        result = zeros( m, ndim);

        v = Nmat{ndim};
        ip = 1;
        for i=1:numel(v),
          isize = index_eq_count( ndim-1,Nmat,Levsum-v(i));
          
          result(ip:(ip+isize-1), 1:(ndim-1)) = ...
                        index_eq( ndim-1,Nmat,Levsum-v(i));
          result(ip:(ip+isize-1),ndim) = i;

          ip = ip + isize;
        end;
end

end
        


