function [flops1, flops2, imethod,imem1,imem2] = flops_kron6( nrow1,ncol1, ...
                                                  nrow2,ncol2, ...
                                                  nrow3,ncol3, ...
                                                  nrow4,ncol4, ...
                                                  nrow5,ncol5, ...
                                                  nrow6,ncol6)
%
% [flops1, flops2, imethod,imem1,imem2] = flops_kron6( nrow1,ncol1, ...
%                                          nrow2,ncol2, ...
%                                          nrow3,ncol3, ...
%                                          nrow4,ncol4, ...
%                                          nrow5,ncol5)
%
% estimate flops (and temporary storage) for method 1 and method 2
%
flops1 = 0;
flops2 = 0;
imem1 = 0;
imem2 = 0;

% -----------------------------------------------
% reshape X as (ncol2*ncol3*ncol4*ncol5*ncol6) by ncol1
% -----------------------------------------------
nrow_X = (ncol2*ncol3*ncol4*ncol5*ncol6);
ncol_X = ncol1;


% ----------------------
% method 1:
% Ytmp = kron(A2,A3,A4,A5,A6)*X
% Y = Ytmp * transpose(A1)
% ----------------------
nvec = ncol_X;
[jflops1,jflops2,jmethod,jmem1,jmem2] = ...
       flops_kron5( nrow2,ncol2,  ...
                    nrow3,ncol3,  ...
                    nrow4,ncol4,  ...
                    nrow5,ncol5,  ...
                    nrow6,ncol6);

if (jmethod == 1),
  flops1 = flops1 + jflops1 * nvec;
  imem1 = imem1 + jmem1 * nvec;
else
  flops1 = flops1 + jflops2 * nvec;
  imem1 = imem1 + jmem2 * nvec;
end;
%
% ------------------------------------------- 
% Ytmp is (nrow2*nrow3*nrow4*nrow5*nrow6) by ncol_X;
% ------------------------------------------- 
nrow_Ytmp = (nrow2*nrow3*nrow4*nrow5*nrow6);  
ncol_Ytmp = ncol_X;

imem1 = imem1 + nrow_Ytmp * ncol_Ytmp;
flops1 = flops1 + (2.0 * nrow_Ytmp ) * ncol1 * nrow1;


% -------------- 
% method 2:
% Ytmp =   X * transpose(A1);
% Y = kron(A2,A3,A4,A5,A6)*Ytmp
% -------------- 

% -------------------------------------
% X reshaped as  nrow_X by ncol_X
%
% Ytmp = X * transpose(A1)
% so  Ytmp is nrow_X * nrow1;
% -------------------------------------
nrow_Ytmp = nrow_X;
ncol_Ytmp = nrow1;

% ------------------------
% Ytmp = X * transpose(A1)
% ------------------------
imem2 = imem2 + nrow_Ytmp * ncol_Ytmp;
flops2 =  flops2 + (2.0 * nrow_X ) * ncol1 * nrow1;


% ----------------------------
% Y = kron(A2,A3,A4,A5) * Ytmp
% ----------------------------
nvec = ncol_Ytmp;
[jflops1,jflops2,jmethod,jmem1,jmem2] = ...
    flops_kron4( nrow2,ncol2,   nrow3,ncol3,   nrow4,ncol4,  nrow5,ncol5, nrow6,ncol6);

if (jmethod == 1),
  flops2 = flops2 +  jflops1 * nvec;
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
