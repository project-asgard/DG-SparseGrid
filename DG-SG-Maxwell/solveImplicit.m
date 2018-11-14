function [xans] = solveImplicit( B, dt, brhs )
% result = solveImplicit( B, dt, brhs )
% solve   
% [ eye   -dt*M ] * xans = brhs
% [ dt*M   eye  ]
%
% where 
% M =  
% [ 0,    M12,  M13;
%  -M12,  0,    M23;
%  -M13,  -M23,  0]
%
% M12 = kron(eye,eye,B)
% M13 = -kron(eye,B,eye)
% M23 = kron(B,eye,eye)
%
idebug = 1;
estimate_spectrum = 1;
n = size(B,1);
isok = (size(B,2) == n);
if (~isok),
   error(sprintf('solveImplicit: B not square,size(B)=%d %d', ...
          size(B,1), size(B,2) ));
  return;
end;

% --------------------------------
% perform eigen decomposition of B
% --------------------------------
if (idebug >= 1),
        disp(sprintf('solveImplicit: size(B)=%d',size(B,1)));
end;
[V,D] = eig(full(B));
% ---------------------------
% create D has sparse diagonal matrix
% ---------------------------
D = sparse(  1:n, 1:n, diag(D), n,n );
E = speye(n,n);
%
% -------------------
% block decomposition
% -------------------
% M12 = kron(eye,eye,B)
% M13 = -kron(eye,B,eye)
% M23 = kron(B,eye,eye)

% ----------------------------------------
% note these are (n^3) by (n^3)  matrices
% ----------------------------------------
D12 =  kron(E,kron(E,D));
D13 = -kron(E,kron(D,E));
D23 =  kron(D,kron(E,E));


% ----------------------------------------------
% Vmat = [V3, 0 ;
%         0,  V3]
%
% V3 = [  kron(V,V,V), 0,           0; 
%         0,           kron(V,V,V), 0;
%         0,           0,           kron(V,V,V)]
% ----------------------------------------------

Z = sparse( n^3,n^3);
Mdiag = [ Z,    D12, D13; ...
         -D12,  Z,   D23; ...
         -D13, -D23, Z];

E3 = speye(3*n^3,3*n^3);
Adiag = [E3,             -dt * Mdiag; ...
         dt * Mdiag,      E3 ];


% ------------------------
% we have
% A = Vmat * Adiag * inv(Vmat) 
%
% solve A * xans = brhs
% Vmat * (Adiag * (inv(Vmat) * xans)) = brhs  
%
% brhs_hat = inv(Vmat) * brhs
% Adiag * yans = brhs_hat,    
% where yans = inv(Vmat)*xans or  xans = Vmat * yans 
% -------------------------

% ---------------------------------------------------
% compute brhs_hat =  inv(Vmat) * brhs
% Note that
% inv( kron(V,V,V) ) = kron( inv(V), inv(V), inv(V) )
% ---------------------------------------------------
isok = (size(brhs,1) == n^3*6);
if (~isok),
   error(sprintf('solveImplicit: incorrect size(brhs)=%d,n=%d,6*n^3=%d', ...
             size(brhs,1), n, 6*n^3 ));
   return;
end;
          
brhs_hat = zeros(size(brhs));
Vinv = inv(V);
for i=1:6,
  i1 = 1 + (i-1)*n^3;
  i2 = i1 + n^3-1;
  brhs_hat(i1:i2) = kron_mult3( Vinv, Vinv, Vinv, brhs(i1:i2));
end;        

% --------------------------------------
% solve  Adiag * yans = brhs_hat
% permute Adiag to be 6x6 block diagonal
% --------------------------------------
iperm = reshape(transpose(reshape(1:(6*n^3), n^3, 6)), 6*n^3,1);
%
% --------------------------------
% P*Adiag*P * P'*yans = P*brhs_hat
% --------------------------------
yans = zeros(6*n^3,1);
use_kxk = 1;
k = 6;
if (use_kxk),
  for iblock=1:(size(Adiag,1)/k),
      i1 = 1 + (iblock-1)*k;
      i2 = i1 + k-1;
      Ai = Adiag(iperm(i1:i2), iperm(i1:i2));
      bi = brhs_hat(iperm(i1:i2));
      yi = Ai\bi;
      yans(iperm(i1:i2)) = yi;
  end;



if (estimate_spectrum),
  nblocks = size(Adiag,1)/k;
  disp(sprintf('nblocks = %d', nblocks));

  lamb = zeros(6,nblocks);
  for iblock=1:(size(Adiag,1)/k),
      i1 = 1 + (iblock-1)*k;
      i2 = i1 + k-1;
      Ai = Adiag(iperm(i1:i2), iperm(i1:i2));
      lamb(1:6,iblock) = eig( full(Ai) );
  end;

  lamb = reshape( lamb, prod(size(lamb)),1);

  figure;
  plot( real( lamb(:)), imag(lamb(:)),'.');
  title('spectrum of A');

  min_re = min( real(lamb));
  max_re = max( real(lamb));
  min_im = min( imag(lamb));
  max_im = max( imag(lamb));
  min_abs = min( abs(lamb));
  max_abs = max( abs(lamb));
  disp(sprintf(...
      'min_re=%g,max_re=%g,min_im=%g,max_im=%g,min_abs=%g,max_abs=%g',...
       min_re,   max_re,   min_im,   max_im,   min_abs,   max_abs));

end;





else
  yans(iperm) = Adiag( iperm, iperm) \ brhs_hat(iperm);
end;
if (idebug >= 2),
        figure;
        spy( Adiag(iperm,iperm) ); 
        title('Adiag(iperm,iperm)');
end;

% ---------------------------
% recover xans as Vmat * yans
% ---------------------------
xans = zeros( size(yans) );
for i=1:6,
   i1 = 1 + (i-1)*n^3;
   i2 = i1 + n^3 - 1;
   xans(i1:i2) = kron_mult3( V,V,V, yans(i1:i2) );
end;


% --------------
% check solution
%
% res = A * xans - brhs
% --------------
if (idebug >= 1),


  % ------------
  % evaluate A*xans
  % ------------
  x1 = xans( 1:(3*n^3) );
  x2 = xans( (3*n^3) + (1:(3*n^3)) );
  Ax1 =  x1              + (-dt)*multM( B,x2);
  Ax2 =  dt*multM(B, x1) + x2;
  Ax = [Ax1; ...
        Ax2];
  resid = brhs - Ax;
  disp(sprintf('norm(resid,1)=%g, norm(brhs,1)=%g, norm(xans,1)=%g', ...
                norm(resid,1),    norm(brhs,1),    norm(xans,1) ));
  
  
  if (idebug >= 2),
    % -------------
    % form A and M
    % -------------
    Z = sparse(n^3,n^3);
    M = formM( B );
    
    eye3 = speye(3*n^3,3*n^3);
    A = [eye3,  -dt*M; ...
        dt*M,   eye3 ];              
    
    Ax = A * xans;
    resid = brhs - Ax;
    disp(sprintf('norm(resid,1)=%g, norm(brhs,1)=%g, norm(xans,1)=%g', ...
                  norm(resid,1),    norm(brhs,1),    norm(xans,1) ));


  end;
  
end;

end
