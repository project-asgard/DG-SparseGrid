function Y = kronmult4(A1,A2,A3,A4, X )
% Y = kronmult4(A1,A2,A3,A4, X )

nrow1 = size(A1,1);
ncol1 = size(A1,2);

nrow2 = size(A2,1);
ncol2 = size(A2,2);

nrow3 = size(A3,1);
ncol3 = size(A3,2);

nrow4 = size(A4,1);
ncol4 = size(A4,2);

nrowX = ncol1*ncol2*ncol3*ncol4;
nvec = numel( X )/nrowX;

nrowY = nrow1*nrow2*nrow3*nrow4;


Ytmp = kronmult3( A2,A3,A4, X );

Ytmp = reshape(Ytmp, [numel(Ytmp)/nvec, nvec]);
nrowYtmp = size(Ytmp,1);

Y = zeros(nrowY, nvec);
for i=1:nvec,
 Yi = reshape( Ytmp(:,i), [nrowYtmp/ncol1, ncol1]) * transpose(A1);
 Y(:,i) = reshape( Yi, [nrowY,1]);
end;

