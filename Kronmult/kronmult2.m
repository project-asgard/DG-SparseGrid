function Y = kronmult2(A1,A2, X )
% Y = kronmult2(A1,A2, X )
global idebug;

nrow1 = size(A1,1);
ncol1 = size(A1,2);

nrow2 = size(A2,1);
ncol2 = size(A2,2);

nrowX = ncol1*ncol2;

isok = (mod(numel(X),nrowX) == 0);
if (~isok),
  error(sprintf('kronmult2: numel(X)=%g, nrowX=%g', ...
                            numel(X),    nrowX ));
  return;
end;

nvec = numel(X)/nrowX;

nrowY = nrow1*nrow2;
Y = zeros( nrowY, nvec );

[flops1,flops2,imethod] = flops_kron2( nrow1,ncol1, nrow2,ncol2);

if (idebug >= 1),
  disp(sprintf('kronmult2: flops1=%g, flops2=%g, imethod=%d', ...
                           flops1,    flops2,    imethod ));
end;


use_method_1 = (imethod == 1);
if (use_method_1),

  Ytmp  = kronmult1( A2, X );
  if (idebug >= 1),
    disp(sprintf('kronmult2: numel(Ytmp)=%g', numel(Ytmp)));
  end;
  
  Ytmp = reshape( Ytmp, [numel(Ytmp)/nvec,nvec]);
  nrowYtmp = size(Ytmp,1);
  
  
  use_single_call = 0;
  if (use_single_call),
    % -------------------------------------------
    % note: just change of view, no data movement
    % -------------------------------------------
    isok = (mod(numel(Ytmp), (ncol1*nvec)) == 0);
    if (~isok),
      error(sprtinf('kronmult2: numel(Ytmp)=%g, ncol1=%g, nvec=%g', ...
                                numel(Ytmp),    ncol1,    nvec ));
      return;
    end;
  
  
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
    isok = (mod(nrowYtmp,ncol1) == 0);
    if (~isok),
     error(sprintf('kronmult2: nrowYtmp=%g, ncol1=%g', ...
                            nrowYtmp,    ncol1 ));
     return;
    end;
  
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
              
  end; % use_single_gemm
else
%   -----------------------------
%   Y = A2 * (X * transpose(A1) )
%   -----------------------------
    X = reshape( X, ncol1*ncol2, nvec );
    Ytmp = zeros( ncol2, nrow1 * nvec );
    for i=1:nvec,
         Xi = reshape( X(:,i), ncol2, ncol1 );
         i1 = 1 + (i-1)*nrow1;
         i2 = i1 + nrow1 - 1;
         
         Ytmp(:,i1:i2) =  Xi(1:ncol2,1:ncol1) * ...
                             transpose( A1(1:nrow1,1:ncol1));
    end;
    Y =  A2(1:nrow2,1:ncol2) * Ytmp(1:ncol2, 1:(nrow1*nvec));
    Y = reshape(  Y, (nrow2*nrow1), nvec );
end;


Y = reshape(  Y, (nrow2*nrow1), nvec );

end;


