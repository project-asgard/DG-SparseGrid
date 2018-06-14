function [flops1, flops2, imethod] = flops_kron5( nrow1,ncol1, ...
                                                  nrow2,ncol2, ...
                                                  nrow3,ncol3, ...
                                                  nrow4,ncol4, ...
                                                  nrow5,ncol5)
%
% [flops1, flops2, imethod] = flops_kron5( nrow1,ncol1, ...
%                                          nrow2,ncol2, ...
%                                          nrow3,ncol3, ...
%                                          nrow4,ncol4, ...
%                                          nrow5,ncol5)
%
% estimate flops for method 1 and method 2
%

% ----------------------
% method 1:
% Ytmp = kron(A2,A3,A4,A5)*X
% Y = Ytmp * transpose(A1)
% ----------------------
nvec = ncol1;
flops1 = flops_kron4( nrow2,ncol2, nrow3,ncol3, nrow4,ncol4,nrow5,ncol5)*nvec;
%
% Ytmp is (nrow2*nrow3*nrow4*nrow5) by ncol1;
% 
nrow_Ytmp = (nrow2*nrow3*nrow4*nrow5);  
flops1 = flops1 + 2.0 * nrow_Ytmp * ncol1 * nrow1;


% -------------- 
% method 2:
% Ytmp =   X * transpose(A1);
% Y = kron(A2,A3)*Ytmp
% -------------- 

% -------------------------------------
% X reshaped as  (ncol2*ncol3*ncol4*ncol5) by ncol1
% Ytmp is (ncol2*ncol3*ncol4*ncol5) by nrow1
% -------------------------------------
flops2 =  2.0 * (ncol2*ncol3*ncol4*ncol5) * ncol1 * nrow1;

flops2 = flops2 + flops_kron4( nrow2,ncol2, nrow3,ncol3, ...
                               nrow4,ncol4, nrow5,ncol5)*nrow1;

if (flops1 <= flops2),
  imethod = 1;
else
  imethod = 2;
end;
