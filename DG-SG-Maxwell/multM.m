function [Mw] = multM( B , w )
% function [Mw] = multM( B , w )
% compute   Mw = M * w
% where
% --------------------
% M = [ 0,    M12,  M13;
%      -M12,  0,    M23;
%      -M13, -M23,  0]
%
% ---------------------
% M12 = kron(eye,eye,B)
% M13 = -kron(eye,B,eye)
% M23 = kron(B,eye,eye)
% ---------------------

idebug = 0;
n = size(B,1);
E = speye(n,n);


% ---------------------
% compute M * w
% M =  
% [ 0,    M12,  M13;
%  -M12,  0,    M23;
%  -M13,  -M23,  0]
%
% M12 =  kron(eye,eye,B)
% M13 = -kron(eye,B,eye)
% M23 =  kron(B,eye,eye)
% ---------------------

ip = 1:(n^3);
w1 = w( ip,:);
w2 = w( n^3 + ip,:);
w3 = w( 2*n^3 + ip, :);

% -----------------------------------
% M12w2 = M12 * w2 =  kron(eye,eye,B)*w2
% M13w3 = M13 * w3 = -kron(eye,B,eye)*w3
% -----------------------------------
M12w2 =  kron_mult3( E,E,B, w2);
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

% ------------
% double check
% ------------
if (idebug >= 1),
  M = formM( B );
  err = norm( M*w - Mw,1);
  disp(sprintf('multM: err=%g,norm(M,1)=%g, norm(w,1)=%g,norm(Mw,1)=%g', ...
                       err,   norm(M,1),    norm(w,1),   norm(Mw,1) ));
end;

