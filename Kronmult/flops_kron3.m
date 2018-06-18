function [flops1, flops2, imethod,imem1,imem2] = flops_kron3( nrow1,ncol1, ...
                                                  nrow2,ncol2, ...
                                                  nrow3,ncol3)

% [flops1, flops2, imethod,imem1,imem2] = flops_kron3( nrow1,ncol1, ...
%                                          nrow2,ncol2, ...
%                                          nrow3,ncol3)
%
% estimate flops (and temporary storage) for method 1 and method 2
%
flops1 = 0;
flops2 = 0;
imem1 = 0;
imem2 = 0;

% -----------------------------------
% reshape X as (ncol2*ncol3) by ncol1;
% -----------------------------------
nrow_X = ncol2*ncol3;
ncol_X = ncol1;

% ----------------------
% method 1:
% Ytmp = kron(A2,A3)*X
% Y = Ytmp * transpose(A1)
% ----------------------


% --------------------
% Ytmp = kron(A2,A3)*X
% --------------------
nvec = ncol_X;
[jflops1, jflops2, jmethod, jmem1, jmem2] = flops_kron2(nrow2,ncol2,nrow3,ncol3);
if (jmethod == 1),
  flops1 = flops1 +  jflops1 * nvec;
  imem1 = imem1 + jmem1 * nvec;
else
  flops1 = flops1 + jflops2 * nvec;
  imem1 = imem1 + jmem2 * nvec;
end;
% ------------------------------- 
% Ytmp is (nrow2*nrow3) by ncol1;
% ------------------------------- 
nrow_Ytmp = (nrow2*nrow3);  
ncol_Ytmp = ncol_X;



% ------------------------
% Y = Ytmp * transpose(A1)
% ------------------------
imem1 = imem1 +  nrow_Ytmp * ncol_Ytmp;
flops1 = flops1 + (2.0 * nrow_Ytmp) * ncol1 * nrow1;


% -------------- 
% method 2:
% Ytmp =   X * transpose(A1);
% Y = kron(A2,A3)*Ytmp
% -------------- 


% --------------------
% Ytmp = X * transpose(A1)
% so Ytmp is nrow_X by nrow1
% --------------------
nrow_Ytmp = nrow_X;
ncol_Ytmp = nrow1;

flops2 =  flops2 + (2.0 * nrow_X) * ncol1 * nrow1;
imem2 = imem2 + nrow_Ytmp * ncol_Ytmp;


% --------------------
% Y = kron(A2,A3)*Ytmp
% --------------------
nvec = ncol_Ytmp;
[jflops1, jflops2, jmethod, jmem1, jmem2 ] = flops_kron2( nrow2,ncol2, nrow3,ncol3);
if (jmethod == 1),
  flops2 = flops2 + jflops1 * nvec;
  imem2 = imem2 + jmem1 * nvec;
else
  flops2 = flops2 + jflops2 * nvec;
  imem2 = imem2 + jmem2 * nvec;
end;


if (flops1 <= flops2),
  imethod = 1;
else
  imethod = 2;
end;
