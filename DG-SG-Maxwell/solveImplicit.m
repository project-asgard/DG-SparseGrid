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
yans(iperm) = Adiag( iperm, iperm) \ brhs_hat(iperm);

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

itop = 1:(3*n^3);
ibot = (3*n^3) + itop;
xtop = xans(itop);
xbot = xans(ibot);

% ---------------------------
% A = [eye ,  -dt*M] * [xtop;
%     [dt*M,   eye ]    xbot]
% ---------------------------
Ax = zeros(size(brhs));

Ax(itop) = xans(itop);
Ax(ibot) = xans(ibot);


w = [xtop, xbot];
% compute M * w
% M =  
% [ 0,    M12,  M13;
%  -M12,  0,    M23;
%  -M13,  -M23,  0]
%
% M12 = kron(eye,eye,B)
% M13 = -kron(eye,B,eye)
% M23 = kron(B,eye,eye)

w1 = w( 1:(n^3),: );
w2 = w( (n^3) + (1:(n^3)),: );
w3 = w( 2*n^3 + (1:(n^3)),: );
% -----------------------------------
% M12w2 = M12 * w2 = kron(eye,eye,B)*w2
% M13w3 = M13 * w3 = -kron(eye,B,eye)*w3
% -----------------------------------
M12w2 =  kron_mult3( E,E,B, w1);
M13w3 = -kron_mult3( E,B,E, w3);
Mw1 = M12w2 + M13w3;
clear M12w2
clear M13w3

% -----------------------------------
% M12w1 = M12*w1 = kron(eye,eye,B)*w1
% M23w3 = M23*w3 = kron(B,eye,eye)*w3
% -----------------------------------
M12w1 = kron_mult3( E,E,B, w1 );
M23w3 = kron_mult3( B,E,E, w3 );
Mw2 = -M12w1 + M23w3;
clear M12w1;
clear M13w3;

% -----------------------------------
% M13w1 = M13*w1 = -kron(eye,B,eye)*w1
% M23w2 = M23*w2 =  kron(B,eye,eye)*w2
% -----------------------------------
M13w1 = -kron_mult3(E,B,E, w1);
M23w2 =  kron_mult3(B,E,E, w2);
Mw3 = -M13w1 - M23w2;
clear M13w1;
clear M23w2;
Mw = [Mw1; ...
      Mw2; ...
      Mw3];

clear Mw1
clear Mw3
clear Mw3

Mxtop = Mw(1:(3*n^3),1);
Mxbot = Mw(1:(3*n^3),2);


Ax(itop) = Ax(itop) - dt*Mxbot;
Ax(ibot) = Ax(ibot) + dt*Mxtop;

resid = brhs - Ax;
disp(sprintf('norm(resid,1)=%g, norm(brsh,1)=%g, norm(xans,1)=%g', ...
              norm(resid,1),    norm(brhs,1),    norm(xans,1) ));



% -------------
% form A and M
% -------------
Z = sparse(n^3,n^3);
% ---------------------
% M12 = kron(eye,eye,B)
% M13 = -kron(eye,B,eye)
% M23 = kron(B,eye,eye)
% ---------------------
M12 = kron(E,kron(E,B));
M13 = -kron(E,kron(B,E));
M23 = kron(B,kron(E,E));

M = [Z,     M12, M13; ...
     -M12,  Z,   M23; ...
     -M13, -M23, Z ];
A = [speye(3*n^3,3*n^3),  -dt*M; ...
    dt*M,                 speye(3*n^3,3*n^3)];

Ax = A * xans;
resid = brhs - Ax;
disp(sprintf('norm(resid,1)=%g, norm(brsh,1)=%g, norm(xans,1)=%g', ...
              norm(resid,1),    norm(brhs,1),    norm(xans,1) ));

end;
