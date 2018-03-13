nvec = 5;
maxit = 100;
nerr = 0;
total_flops_performed = 0;
idebug = 0;
total_time = 0;

for it=1:maxit,
  idim = round( rand(4,1) * 100 ) + 1;
  nrowA = idim(1); ncolA = idim(2);
  nrowB = idim(3); ncolB = idim(4);
  A = rand(nrowA,ncolA);
  B = rand(nrowB,ncolB);
  X = rand(ncolB*ncolA, nvec);

  Y = kron(A,B)*X;

  t1 = tic;
  [Y2,flops_performed] = kron_mult2(A,B,X,idebug);
  total_time = total_time + toc(t1);

  total_flops_performed = total_flops_performed + flops_performed;

  Y = reshape( Y, nrowB*nrowA,nvec);
  Y2 = reshape( Y2, nrowB*nrowA,nvec);

  diff = norm(Y-Y2,1);
  normA = norm(A,1);
  normB = norm(B,1);
  isok = abs(diff < 1e-6*max(normA,normB));
  if (~isok),
    disp(sprintf('it=%g, diff=%g, normA=%g, normB=%g', ...
                  it,    diff,    normA,    normB ));
    disp(sprintf('size(A)=(%d,%d), size(B)=(%d,%d) ', ...
                  size(A,1),size(A,2), size(B,1),size(B,2) ));
    nerr = nerr + 1;
   end;
end;
 
if (nerr == 0),
  disp(sprintf('ALL OK after %g tests', maxit)); 
  gflops = total_flops_performed * 10^(-9)/total_time;

  disp(sprintf('overall rate is effectively %g Gflops/sec', gflops ));
else
  disp(sprintf('%d errors ', nerr ));
end;
