function [Y,flops_performed] = kron_mult2(A,B,X, idebug_in)
% [Y, flops_performed]  = kron_mult2(A,B,X)
%
idebug = 0;
if (nargin >= 4),
  idebug = idebug_in;
end;

flops_performed = 0;

use_kron_multd = 0;
if (use_kron_multd),
  Acell{1} = A;
  Acell{2} = B;
  nkron = 2;
  Y = kron_multd(nkron,Acell,X);
else
  % -----------------------------------
  % Y = B * reshape(X) * transpose(A)
  % -----------------------------------
  nrowA = size(A,1); ncolA = size(A,2);
  nrowB = size(B,1); ncolB = size(B,2);
  nrowX = ncolB; ncolX = ncolA;
  nrowY = nrowB; ncolY = nrowA;

  isok = mod( numel(X), nrowX*ncolX ) == 0;
  if (~isok),
    error(sprintf('kron_mult2: invalid size of X, numel(X)=%g, nrowX=%g, ncolX=%g', ...
                   numel(X), nrowX, ncolX ));
    return;
  end;
  nvec = numel(X)/(nrowX*ncolX);
  
  nnzA = nnz(A);
  nnzB = nnz(B);


  if (idebug >= 1),
    disp(sprintf('nvec=%d, size(A)=(%d,%d), nnzA=%g, size(B)=(%d,%d),nnzB=%g', ...
                  nvec, ...
                  size(A,1),size(A,2),nnzA, ...
                  size(B,1),size(B,2),nnzB ));
  end;
  
  % -----------------------------------
  % method 1:   
  % BX = B * X;   
  % Y = BX * transpose(A)
  % -----------------------------------
  
  nrowBX = nrowB; ncolBX = ncolX;
  flops_BX = 2.0*nnzB*ncolBX;

  flops_BX_At = 2.0*nrowBX * nnzA;

  flops_method(1) = flops_BX + flops_BX_At;

  % -----------------------------------
  % method 2:
  % XAt = X * transpose(A);
  % Y = B * XAt
  % -----------------------------------
  ncolXAt = nrowA;
  nrowXAt = nrowX;

  flops_XAt = 2.0*nrowXAt * nnzA;

  flops_B_XAt = 2.0*nnzB * ncolXAt;

  flops_method(2) = flops_XAt  + flops_B_XAt;

  % ------------------
  % method 3:
  % expand kron(A,B)*X
  % ------------------
  flops_method(3) = 2.0*nnzA*nnzB;

  min_flops = min( flops_method(1), min( flops_method(2), flops_method(3) ) );
  if (idebug >= 1),
    disp(sprintf('flops_method(:)=(%g,%g,%g) ', ...
           flops_method(1), flops_method(2), flops_method(3)));
  end;


  nvec = numel(X)/(nrowX*ncolX);
  X = reshape( X, [ nrowX, ncolX, nvec ] );
  Y = zeros( [nrowY, ncolY, nvec ] );
  flops_performed = min_flops * nvec;

  if (flops_method(3) <=  min_flops),

    Y = kron(A,B)*reshape( X, nrowX*ncolX, nvec );

  elseif (flops_method(2) <=  min_flops),

     for k=1:nvec,
        XAt = X(1:nrowX,1:ncolX,k) * transpose(A);
        Y(1:nrowY,1:ncolY,k) = B * XAt;
      end;
     Y = reshape( Y, nrowY*ncolY,nvec);

  else

     for k=1:nvec,
         BX = B * X(1:nrowX,1:ncolX,k);
         Y(1:nrowY,1:ncolY,k) = BX* transpose(A);
     end;
     Y = reshape( Y, [nrowY*ncolY,nvec]);
  end;

  

end;
