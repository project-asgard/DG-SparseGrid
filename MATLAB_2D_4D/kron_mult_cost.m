function [flops, schedule] = kron_mult_cost( rc, nz, xsizes )
% [flops, schedule] = kron_mult_cost( rc, nz, xsizes )
% estimate number of flops to perform kronecker multiplication
%  kron(A1,A2,..,Ak) * X, 
% where  size(Aj) is   rc(1,j) by rc(2,j) and has nz(j) non-zeros
%
idebug = 0;
nkron = length(nz);
isok = (size(rc,2) == nkron);
if (~isok),
  error(sprintf('kron_mult_cost: mismatch in sizes of rc (%g,%g) and nz (%g) ',...
               size(rc,1), size(rc,2), length(nz) ));
  return;
end;

isok = (mod( prod(xsizes), prod(rc(2,:)) ) == 0);
if (~isok),
  error(sprintf('kron_mult_cost: mismatch in xsizes %g  and column sizes %g', ...
           prod(xsizes),  prod(rc(2,:))  ));
  return;
end;

if (idebug >= 1),
  disp(sprintf('kron_mult_cost:nkron=%g, prod(xsizes)=%g ', ...
                nkron, prod(xsizes) ));
  for k=1:nkron,
    disp(sprintf('k=%g, rc(1,k)=%g, rc(2,k)=%g ', ...
                  k,    rc(1,k),    rc(2,k) ));
  end;
end;



if (nkron == 1),
  % -----------
  % just A1 * X
  % -----------
  ncolA = rc(2,1);
  nvec = prod(xsizes)/ncolA;
  flops = 2.0 * nz(1) * nvec;
  schedule = [1];

  if (idebug >= 1),
    disp(sprintf('kron_mult_cost: nkron=1, flops=%g ', flops));
  end;
elseif (nkron == 2),
  % ----------------------------------
  %  kron(A1,A2) * X =>  A2*X*A1^t
  %  method 1:
  %
  %    A2X = A2*X  need 2*nnz(A2)*ncolA1 flops
  %    A2X is  nrowA2 by ncolA1
  %
  %    A2X * A1^t or   { A1*(A2X^t) }^t need 2*nnz(A1)*nrowA2
  %
  % ----------------------------------
  nrowA1 = rc(1,1);
  ncolA1 = rc(2,1);
  nrowA2 = rc(1,2);
  ncolA2 = rc(2,2);

  ncolX = prod(xsizes)/(ncolA2*ncolA1);
  
  nnz_A1 = nz(1);
  nnz_A2 = nz(2);
  flops_method1 = 2.0*nnz_A2*ncolA1 + ...
                  2.0*nnz_A1*nrowA2;

  % ----------------------------------
  %   XA1t = X*A1^t  or {A1*X^t}^t
  %   need  2*nnz(A1) * ncolA2
  %
  %   A2 * XA1t need  2*nnz(A2)*ncolA1
  % ----------------------------------
  flops_method2 = 2.0*nnz_A1 * ncolA2 + ...
                  2.0*nnz_A2 * ncolA1;

  % -----------------------------------------
  % expand kron(A1,A2), then perform multiply
  % -----------------------------------------
  flops_method3 = 2.0 * nnz_A1*nnz_A2;

  flops = min(flops_method1,min(flops_method2,flops_method3))*ncolX;
  schedule = [1];

  if (idebug >= 1),
    disp(sprintf(...
     'kron_mult_cost: flops_method1=%g, flops_method2=%g, flops=%g', ...
      flops_method1, flops_method2, flops ));
 end;
             
elseif (nkron >= 3),
  flopsvec = zeros(nkron-1,1);
  
  for k=1:(nkron-1),
    [flops1,schedule1] = kron_mult_cost( rc(:,1:k), nz(1:k), xsizes ); 
    [flops2,schedule2] = kron_mult_cost( rc(:, (k+1):nkron), nz( (k+1):nkron), xsizes );
    sch{k} = [schedule1,schedule2];
    flopsvec(k) = flops1 + flops2;
  end;

  if (idebug >= 1),
   for k=1:(nkron-1),
     disp(sprintf('kron_mult_cost:k=%g, flops=%g', ...
               k, flopsvec(k) ));
   end;
  end;

  [min_flops, min_k] = min( flopsvec );
  flops = min_flops;

  % -------------------------------------
  % consider case of very sparse matrices
  % -------------------------------------
  sparse_flops = 2.0*prod(nz(1:nkron));
  if (sparse_flops < min_flops),
    schedule = [0];
  else 
    schedule = [sch{min_k}];
  end;
end;

  
   
                  


