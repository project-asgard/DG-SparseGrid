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

nrowY = nrow1*nrow2*nrow3*nrow4;


% --------------------------------
% Ytmp = kronmult3( A2,A3,A4, X );
% --------------------------------
[batch_list1, batch_list2, batch_list3, Ytmp] = ...
          kronmult3_batch(A2,A3,A4, X,  ...
                  batch_list1, batch_list2, batch_list3 );
if (idebug >= 1),
  disp(sprintf('kronmult4_batch: numel(Ytmp)=%g', numel(Ytmp)));
end;

Ytmp = reshape(Ytmp, [numel(Ytmp)/nvec, nvec]);
nrowYtmp = size(Ytmp,1);

Y = zeros(nrowY, nvec);

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

 Amat = reshape( Ytmp(:,i), mm,kk);
 Bmat = A1;
 Cmat = reshape( Y(:,i), mm, nn );

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

  Y(:,i) = reshape( Cmat, [nrowY,1]);
end;

end
