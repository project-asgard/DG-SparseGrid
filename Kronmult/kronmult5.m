function Y = kronmult5(A1,A2,A3,A4,A5, X)
% Y = kronmult5(A1,A2,A3,A4,A5, X)

nrow1 = size(A1,1);
ncol1 = size(A1,2);

nrow2 = size(A2,1);
ncol2 = size(A2,2);

nrow3 = size(A3,1);
ncol3 = size(A3,2);

nrow4 = size(A4,1);
ncol4 = size(A4,2);

nrow5 = size(A5,1);
ncol5 = size(A5,2);

nrowX = ncol1*ncol2*ncol3*ncol4*ncol5;

isok = (mod(numel(X), nrowX) == 0);
if (~isok),
 error(sprintf('kronmult5: numel(X)=%g, nrowX=%g', ...
                           numel(X),    nrowX ));
 return;
end;


nvec = numel( X )/nrowX;

nrowY = nrow1*nrow2*nrow3*nrow4*nrow5;


Ytmp = kronmult4( A2,A3,A4, A5, X );

isok = (mod(numel(Ytmp),nvec) == 0);
if (~isok),
  error(sprintf('kronmult5: numel(Ytmp)=%g, nvec=%g', ...
                            numel(Ytmp),    nvec ));
  return;
end;

nrowYtmp = numel(Ytmp)/nvec;
Ytmp = reshape(Ytmp, [numel(Ytmp)/nvec, nvec]);

Y = zeros(nrowY, nvec);

for i=1:nvec,
 Yi = reshape( Ytmp(:,i), [nrowYtmp/ncol1, ncol1]) * transpose(A1);
 Y(:,i) = reshape( Yi, [nrowY,1]);
end;

