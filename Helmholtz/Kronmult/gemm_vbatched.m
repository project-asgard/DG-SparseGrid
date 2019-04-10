function [Clist] = gemm_vbatched(transA, transB, mlist, nlist, klist, ...
                                 alpha, Alist, Blist, beta, Clist )
%
% [Clist] = gemm_vbatched(transA, transB, mlist, nlist, klist, ...
%                                 alpha, Alist, Blist, beta, Clist )
%
% note assume all cases have same scalar values for 
% transA, transB, alpha, beta
% to match magmablas_dgemm_vbatched

nbatch = numel(mlist);
isok = (numel(nlist) == nbatch) && ...
       (numel(klist) == nbatch) && ...
       (numel(Alist) == nbatch) && ...
       (numel(Blist) == nbatch) && ...
       (numel(Clist) == nbatch);
if (~isok),
  error(sprintf('gemm_vbatched: mismatch nbatch=%d',nbatch));
  return;
end;

for ibatch=1:nbatch,
   mm = mlist(ibatch);
   nn = nlist(ibatch);
   kk = nlist(ibatch);

   Amat = Alist{ibatch};
   Bmat = Blist{ibatch};
   Cmat = Clist{ibatch};

   Cmat = gemm( transA, transB, mm,nn,kk, ...
                  alpha, Amat,  Bmat,       ...
                  beta,  Cmat );

   Clist{ibatch} = Cmat;
end;

end
