function Y = kronmult3(A1,A2,A3, X )
% Y = kronmult3(A1,A2,A3, X )
global idebug;

nrow1 = size(A1,1);
ncol1 = size(A1,2);

nrow2 = size(A2,1);
ncol2 = size(A2,2);

nrow3 = size(A3,1);
ncol3 = size(A3,2);

nrowX = ncol1*ncol2*ncol3;

isok = (mod( numel(X), nrowX ) == 0);
if (~isok),
 error(sprintf('kronmult3: numel(X)=%g, nrowX=%g', ...
                                numel(X),    nrowX));
 return;
end;

nvec = numel(X)/nrowX;

nrowY = nrow1*nrow2*nrow3;

Ytmp = kronmult2( A2,A3,X);
if (idebug >= 1),
  disp(sprintf('kronmult3: numel(Ytmp)=%g', numel(Ytmp)));
end;

isok = (mod( numel(Ytmp), nvec ) == 0);
if (~isok),
  error(sprintf('kronmult3: numel(Ytmp)=%g, nvec=%g', ...
                            numel(Ytmp),    nvec));
  return;
end;

nrowYtmp = numel(Ytmp)/nvec;
Ytmp = reshape( Ytmp, [nrowYtmp, nvec] );

Y = zeros(nrowY, nvec);

% -----------------------------------------------------
% note: may be task parallelism or batch gemm operation
% -----------------------------------------------------
isok = (mod(nrowYtmp,ncol1) == 0);
if (~isok),
  error(sprintf('kronmult3: nrowYtmp=%g, ncol1=%g', ...
                            nrowYtmp,    ncol1 ));
  return;
end;
msize  = nrowYtmp/ncol1;
for i=1:nvec,
  Yi = reshape(Ytmp(:,i), [msize, ncol1])*transpose(A1);
  Y(:,i) = reshape(Yi, nrowY,1);
end;
