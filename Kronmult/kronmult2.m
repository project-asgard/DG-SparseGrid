function Y = kronmult2(A1,A2, X )
% Y = kronmult2(A1,A2, X )
idebug = 1;

nrow1 = size(A1,1);
ncol1 = size(A1,2);

nrow2 = size(A2,1);
ncol2 = size(A2,2);

nrowX = ncol1*ncol2;
nvec = numel(X)/nrowX;

nrowY = nrow1*nrow2;

Ytmp  = kronmult1( A2, X );
Ytmp = reshape( Ytmp, [numel(Ytmp)/nvec,nvec]);
nrowYtmp = size(Ytmp,1);

Y = zeros( nrowY, nvec );

use_single_call = 0;
if (use_single_call),
  % -------------------------------------------
  % note: just change of view, no data movement
  % -------------------------------------------
  n1 = numel(Ytmp)/(ncol1*nvec);
  Ytmp =  reshape( Ytmp, [n1,ncol1,nvec]);

  
  % ------------------------------
  % need axis permutation
  % ------------------------------
  Yin = permute( Ytmp, [1,3,2]);
  Yin = reshape( Yin, [n1*nvec, ncol1]);

  Yout = Yin * transpose(A1);

  if (idebug >= 1),
    mm = size(Yout,1);
    nn = size(Yout,2);
    kk = ncol1;

    disp(sprintf('kronmult2: single call to gemm of (m,n,k)=(%d,%d,%d)', ...
          mm,nn,kk ));
  end;
  Yout = reshape( Yout, [n1, nvec, nrow1]);

  % ------------------------------
  % need axis permutation
  % ------------------------------
  Y = permute( Yout, [1,3,2]);
  Y = reshape( Y, [nrowY, nvec] );

else
  n1 = nrowYtmp/ncol1;

% ----------------------------------
% note may be batched gemm operation
% ----------------------------------
  for i=1:nvec,
   Yi = reshape( Ytmp(:,i), [n1, ncol1])*transpose(A1);
   Y(:,i) = reshape(Yi, nrowY,1);
  end;

  if (idebug >= 1),
    mm = n1; nn = nrow1; kk = ncol1;
    disp(sprintf('kronmult2: nvec=%d calls to gemm of (m,n,k)=(%d,%d,%d)', ...
                             nvec, mm,nn,kk ));
  end;
            
end;


