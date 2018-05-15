function Y = kronmult3(A1,A2,A3, X )
% Y = kronmult3(A1,A2,A3, X )

nrow1 = size(A1,1);
ncol1 = size(A1,2);

nrow2 = size(A2,1);
ncol2 = size(A2,2);

nrow3 = size(A3,1);
ncol3 = size(A3,2);

nrowX = ncol1*ncol2*ncol3;
nvec = numel(X)/nrowX;

nrowY = nrow1*nrow2*nrow3;

Ytmp = kronmult2( A2,A3,X);
Ytmp = reshape( Ytmp, [numel(Ytmp)/nvec, nvec] );
nrowYtmp = size(Ytmp,1);

Y = zeros(nrowY, nvec);

for i=1:nvec,
  Yi = reshape(Ytmp(:,i), [nrowYtmp/ncol1, ncol1])*transpose(A1);
  Y(:,i) = reshape(Yi, nrowY,1);
end;
