function [batch_list1, batch_list2,Y] = kronmult2_batch(A1,A2, X, batch_list1, batch_list2 )
% [batch_list1, batch_list2,Y] = kronmult2(A1,A2, X, batch_list1, batch_list2 )
% 
global idebug;

nrow1 = size(A1,1);
ncol1 = size(A1,2);

nrow2 = size(A2,1);
ncol2 = size(A2,2);

nrowX = ncol1*ncol2;

isok = (mod(numel(X),nrowX) == 0);
if (~isok),
  error(sprintf('kronmult2_batch: numel(X)=%g, nrowX=%g', ...
                            numel(X),    nrowX ));
  return;
end;

nvec = numel(X)/nrowX;

if (idebug >= 1),
  disp(sprintf('kronmult2_batch: (%d,%d) (%d,%d) nvec=%d', ...
            nrow1,   ncol1,    nrow2,   ncol2,    nvec));
end;

nrowY = nrow1*nrow2;

Y = zeros( nrowY, nvec );

[flops1,flops2,imethod] = flops_kron2( nrow1,ncol1, nrow2,ncol2);
use_method_1  = (imethod == 1);
if (use_method_1),

  % --------------------------
  % Ytmp  = kronmult1( A2, X );
  % --------------------------
  [batch_list1, Ytmp] = kronmult1_batch( A2, X, batch_list1);
  
  % ---------------------
  % Ytmp is (nrow2*ncol1)  by nvec
  % ---------------------
  nrowYtmp = nrow2 * ncol1;
  isok = (numel(Ytmp) == (nrowYtmp * nvec ));
  if (~isok),
    disp(sprintf('kronmult2_batch: numel(Ytmp)=%d,nrowYtmp=%d, nrow2=%d', ...
                                   numel(Ytmp),   nrowYtmp,    nrow2 ));
    disp(sprintf('nrow1=%d,ncol1=%d, nrow2=%d,ncol2=%d', ...
                  nrow1,   ncol1,    nrow2,   ncol2));
    disp(sprintf('numel(X)=%d, nrowX=%d, nvec=%d', ...
                  numel(X),    nrowX,    nvec ));
  end;
  Ytmp = reshape( Ytmp, nrowYtmp, nvec );




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
   % ----------------------------------------------------
   % Yi = reshape( Ytmp(:,i), [n1, ncol1])*transpose(A1);
   % ----------------------------------------------------
   transA = 'N';
   transB = 'T';
   mm = n1;
   kk = ncol1;
   nn = nrow1;
   alpha = 1;  
   beta = 0;
   Amat = reshape( Ytmp(1:(mm*kk),i), mm,kk);
   Bmat = A1;
   Cmat = reshape( Y(1:(mm*nn),i), mm,nn );

   nbatch = batch_list2.nbatch;
   nbatch = nbatch + 1;
   batch_list2.nbatch = nbatch;

   batch_list2.mlist(nbatch) = mm;
   batch_list2.nlist(nbatch) = nn;
   batch_list2.klist(nbatch) = kk;

   batch_list2.transA(nbatch) = transA;
   batch_list2.transB(nbatch) = transB;

   batch_list2.alpha(nbatch) = alpha;
   batch_list2.beta(nbatch) = beta;

   batch_list2.Alist{nbatch} = Amat;
   batch_list2.Blist{nbatch} = Bmat;
   batch_list2.Clist{nbatch} = Cmat;


   % ----------------------------------------
   % perform computation to check correctness
   % ----------------------------------------
   Cmat = gemm( transA, transB, mm,nn,kk, alpha, Amat, Bmat, beta, Cmat );
   
   Y(1:nrowY,i) = reshape(Cmat(1:mm,1:nn), nrowY,1);
  end;

  if (idebug >= 1),
    mm = n1; nn = nrow1; kk = ncol1;
    disp(sprintf('kronmult2_batch: nvec=%d calls to gemm of (m,n,k)=(%d,%d,%d)', ...
                             nvec, mm,nn,kk ));
  end;
            

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
         
%          -------------------------------------
%          Ytmp(:,i1:i2) =  Xi(1:ncol2,1:ncol1) * ...
%                             transpose( A1(1:nrow1,1:ncol1));
%          -------------------------------------
          transA = 'N';
          transB = 'T';
          mm = ncol2;
          nn = nrow1;
          kk = ncol1;
          alpha = 1;
          beta = 0;

          Amat = reshape( X(:,i), ncol2, ncol1 );
          Bmat = A1;
          Cmat = Ytmp(1:ncol2, i1:i2);


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

         Cmat = gemm( transA, transB, mm,nn,kk, alpha, Amat, Bmat, beta, Cmat );
         
         Ytmp(1:mm, i1:i2) = Cmat(1:mm,1:nn);

    end;
%     ---------------------------------------------------------
%     Y =  A2(1:nrow2,1:ncol2) * Ytmp(1:ncol2, 1:(nrow1*nvec));
%     ---------------------------------------------------------
   mm = nrow2;
   nn = nrow1*nvec;
   kk = ncol2;
   transA = 'N';
   transB = 'N';
   alpha = 1;
   beta = 0;
   Amat = A2(1:mm,1:kk);
   Bmat = Ytmp(1:kk,1:nn);
   Cmat = zeros(mm,nn);

   Cmat = gemm( transA, transB, mm,nn,kk, alpha, Amat, Bmat, beta, Cmat );
   Y = Cmat;

   nbatch = batch_list2.nbatch + 1;
   batch_list2.nbatch = nbatch;


   batch_list2.mlist(nbatch) = mm;
   batch_list2.nlist(nbatch) = nn;
   batch_list2.klist(nbatch) = kk;

   batch_list2.transA(nbatch) = transA;
   batch_list2.transB(nbatch) = transB;

   batch_list2.alpha(nbatch) = alpha;
   batch_list2.beta(nbatch) = beta;

   batch_list2.Alist{nbatch} = Amat;
   batch_list2.Blist{nbatch} = Bmat;
   batch_list2.Clist{nbatch} = Cmat;

    Y = reshape(  Y, (nrow2*nrow1), nvec );

end;




end


