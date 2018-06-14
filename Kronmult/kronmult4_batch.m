function [ batch_list1, batch_list2, batch_list3, batch_list4, Y] = ...
              kronmult4_batch(A1,A2,A3,A4, X, ...
                  batch_list1, batch_list2, batch_list3, batch_list4 )
% [batch_list1, batch_list2, batch_list3, batch_list4, Y] = ...
%               kronmult4_batch(A1,A2,A3,A4, X, ...
%                   batch_list1, batch_list2, batch_list3, batch_list4 )
global idebug;

nrow1 = size(A1,1);
ncol1 = size(A1,2);

nrow2 = size(A2,1);
ncol2 = size(A2,2);

nrow3 = size(A3,1);
ncol3 = size(A3,2);

nrow4 = size(A4,1);
ncol4 = size(A4,2);

nrowX = ncol1*ncol2*ncol3*ncol4;

isok = (mod(numel(X), nrowX) == 0);
if (~isok),
  error(sprintf('kronmult4_batch: numel(X)=%g, nrowX=%g', ...
                            numel(X),    nrowX ));
  return;
end;

nvec = numel( X )/nrowX;

if (idebug >= 1),
 disp(sprintf('kronmult4_batch: (%d,%d)  (%d,%d) (%d,%d) (%d,%d) nvec=%d', ...
        nrow1,ncol1, ...
        nrow2,ncol2, ...
        nrow3,ncol3, ...
        nrow4,ncol4, ...
        nvec ));
end;

nrowY = nrow1*nrow2*nrow3*nrow4;
Y = zeros(nrowY, nvec);

[flops1,flops2,imethod] = flops_kron4( nrow1,ncol1, nrow2,ncol2, ...
                                       nrow3,ncol3, nrow4,ncol4);
use_method_1  = (imethod == 1);
if (use_method_1),
  
  % --------------------------------
  % Ytmp = kronmult3( A2,A3,A4, X );
  % --------------------------------
  [batch_list1, batch_list2, batch_list3, Ytmp] = ...
            kronmult3_batch(A2,A3,A4, X,  ...
                    batch_list1, batch_list2, batch_list3 );
  
  % ------------------------------------------
  % X is  (ncol4*ncol3*ncol2*ncol1) by nvec
  %
  % Ytmp = kronmult3( A2,A3,A4,   X), so X appears to be (ncol4*ncol3*ncol2) by (ncol1*nvec)
  % Ytmp is  (nrow4*nrow3*nrow2) by (ncol1*nvec)
  % so Ytmp can be reshaped to be  (nrow4*nrow3*nrow2 * ncol1) by nvec
  % ------------------------------------------
  
  nrowYtmp = (nrow4*nrow3*nrow2)*ncol1;
  Ytmp = reshape(Ytmp, nrowYtmp, nvec );
  
  
  % ----------------------------------------------
  % note: task parallelism or batch gemm operation
  % ----------------------------------------------
  isok = (mod(nrowYtmp,ncol1) == 0);
  if (~isok),
    error(sprintf('kronmult4_batch: nrowYtmp=%g, ncol1=%g', ...
                              nrowYtmp,    ncol1 ));
    return;
  end;
  
  msize = nrowYtmp/ncol1;
  for i=1:nvec,
   % --------------------------------------------------------
   % Yi = reshape( Ytmp(:,i), [msize, ncol1]) * transpose(A1);
   % --------------------------------------------------------
   mm = msize;
   kk = ncol1;
   nn = nrow1;
   transA = 'N';
   transB = 'T';
   alpha = 1;
   beta = 0;
  
   Amat = reshape( Ytmp(1:(mm*kk),i), mm,kk);
   Bmat = A1;
   Cmat = reshape( Y(1:(mm*nn),i), mm, nn );
  
    nbatch =  batch_list4.nbatch;
    nbatch = nbatch+1;
    batch_list4.nbatch = nbatch;
    
    batch_list4.mlist(nbatch) = mm;
    batch_list4.nlist(nbatch) = nn;
    batch_list4.klist(nbatch) = kk;
  
    batch_list4.transA(nbatch) = transA;
    batch_list4.transB(nbatch) = transB;
  
    batch_list4.alpha(nbatch) = alpha;
    batch_list4.beta(nbatch) = beta;
  
    batch_list4.Alist{nbatch} = Amat;
    batch_list4.Blist{nbatch} = Bmat;
    batch_list4.Clist{nbatch} = Cmat;
  
  % ----------------------------------------
  % perform computation to check correctness
  % ----------------------------------------
    Cmat = gemm( transA, transB, mm, nn,kk, alpha, Amat, Bmat, beta, Cmat);
  
    Y(1:nrowY,i) = reshape( Cmat(1:mm,1:nn), [nrowY,1]);
  end; % end for

else

% -------------------------------------------------
% Y = kron( A2, A3, A4) * (  X * transpose( A1 ) )
% -------------------------------------------------
  X = reshape(X,  nrowX, nvec );
  Ytmp = zeros( ncol2*ncol3*ncol4, nvec * nrow1 );

  for i=1:nvec,
     i1 = 1 + (i-1)*nrow1;
     i2 = i1 + nrow1-1;

%      ------------------------------------------
%      Xi = reshape( X(:,i), ncol2*ncol3*ncol4, ncol1 );
%      Ytmpi(1:(ncol2*ncol3*ncol4), 1:nrow1 ) = ...
%            Xi(1:(ncol2*ncol3*ncol4), 1:ncol1) * ...
%                 transpose( A1(1:nrow1,1:ncol1));
%      Ytmp(1:(ncol2*ncol3*ncol4),i1:i2) = Ytmpi(1:(ncol2*ncol3*ncol4),1:nrow1);
%      ------------------------------------------

       transA = 'N';
       transB = 'T';
       mm = (ncol2*ncol3*ncol4);
       nn = nrow1;
       kk = ncol1;
       alpha = 1;
       beta = 0;
       Amat = reshape( X(:,i), (ncol2*ncol3*ncol4), ncol1 );
       Bmat = A1;
       Cmat = Ytmp(1:mm, i1:i2 );

       Cmat = gemm( transA, transB, mm,nn,kk, alpha, Amat, Bmat, beta, Cmat );
       Ytmp(1:mm, i1:i2) = Cmat(1:mm,1:nn);


       nbatch = batch_list1.nbatch + 1;
       batch_list1.nbatch = nbatch;
    
    
       batch_list1.mlist(nbatch) = mm;
       batch_list1.nlist(nbatch) = nn;
       batch_list1.klist(nbatch) = kk;
    
       batch_list1.transA(nbatch) = transA;
       batch_list1.transB(nbatch) = transB;
    
       batch_list1.alpha(nbatch) = alpha;
       batch_list1.beta(nbatch) = beta;
    
       batch_list1.Alist{nbatch} = Amat;
       batch_list1.Blist{nbatch} = Bmat;
       batch_list1.Clist{nbatch} = Cmat;

  end; % end for

%   --------------------------------
%   Y = kronmult3( A2,A3,A4,  Ytmp );
%   --------------------------------
[batch_list2,batch_list3,batch_list4, Y] = ...
        kronmult3_batch( A2,A3,A4,  Ytmp,  ...
                  batch_list2, batch_list3, batch_list4 );
end

Y = reshape( Y, nrowY, nvec );

end
