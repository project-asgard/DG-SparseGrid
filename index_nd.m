function result = index_nd(Dim, Kstart,Kend, sizes )
%
% result = index_nd(Dim, Kstart,Kend, sizes )
% similar to sub2ind, return linear index
% from n-dimensional subscripts
%
% Note
% Kstart(1:Dim), Kend(1:Dim), sizes(1:Dim)
% --------------------------------------------
if (Dim == 1),
        K_start = Kstart(1);
        K_end = Kend(1);
        result = reshape( K_start:K_end, K_end-K_start+1,1);
else
        Dimm1 = Dim-1;
        res = index_nd( Dimm1, Kstart(1:Dimm1),Kend(1:Dimm1), sizes(1:Dimm1));
        size_res = prod(size(res));
        res = reshape( res, size_res,1);

        K_start = Kstart(Dim);
        K_end = Kend(Dim);
        K_size = K_end - K_start + 1;
        result = zeros( size_res * K_size,1);

        lsize = prod(sizes(1:(Dim-1)));
        for i=1:K_size,
           i1 = (i-1)*size_res + 1;
           i2 = i1 + size_res-1;
           K_value = K_start + (i-1);
           result( i1:i2 ) = res + (K_value - 1) * lsize;
        end;
end;
return;
end
