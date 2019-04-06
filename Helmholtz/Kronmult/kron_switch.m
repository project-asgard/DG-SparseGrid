function [ip_r,ip_c] = kron_switch( nrA,ncA, nrB,ncB )
% [ip_r, ip_c] = kron_switch( nrA,ncA, nrB,ncB )
%
% -----------------------------------------
% Pr * kron(A,B) * Pc = kron(B,A)
%
% C = kron(A,B),  then C( [ib,ia], [jb,ja])
% D = kron(B,A),  then D( [ia,ib], [ja,jb])
% -----------------------------------------

ip_r = transpose(reshape( 1:(nrB*nrA), nrB,nrA));
ip_r = reshape( ip_r, numel(ip_r),1);

ip_c = transpose(reshape( 1:(ncB*ncA), ncB,ncA));
ip_c = reshape( ip_c, numel(ip_c),1);

return
end


