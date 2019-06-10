function FMWT = OperatorTwoScale_nonwavelet(pde,maxDeg,maxLev)
% FMWT = OperatorTwoScale_nonwavelet(maxDeg,maxLev)
%----------------------------------
% Set-up Two-scale operator       %
%----------------------------------
% Input: Degree: maxDeg
%        Level: maxLev
% Output: Convert Matrix: FMWT
%**********************************
% http://amath.colorado.edu/~beylkin/papers/A-B-C-R-1993.pdf
%
k = maxDeg;
L = round(log2(maxLev));
isok = (2^L == maxLev);
if (~isok),
    error('OperatorTwoScale_nonwavelet: maxLev=%d must be a power of 2', ...
                                    maxLev);
    return;
end;

n = k * 2^L;

need_shift_x = 1;
% --------------------------------------------------
% note block{0} or block{-1} is invalid, since
% indexing must start from 1, so provide an offset
% so that  block{ioff + (-1)} is valid
% --------------------------------------------------
ioff = 2;

% ---------------------------------------
% following the paper in A-B-C-R-1993.pdf
% Given  x(1) ... x(n)
% ---------------------------------------
xmin = 0;
xmax = 1;
x = linspace(xmin,xmax,n);

for lev=(L-1):-1:0,
   ncells = 2^lev;
   isize = n/ncells;

   % ------------------------------
   % map [x(1)..x(isize)] to [-1,1]
   % ------------------------------
   xx = x(1:isize);
   xx = (xx - min(xx))/( max(xx)-min(xx) );
   if (need_shift_x),
     xx = 2*xx-1;
   end;

   ndeg = 2*k-1;
   F = legendrepoly( ndeg, xx);
   % ----------------------------
   % apply orthogonalization from 
   % previous levels
   %
   % Projection operation is 
   % Y = (eye - Q*Q')*X
   % or  Y = X - Q * (Q'*X)
   % then Q'*Y = 0
   %
   % check Q'*Y = Q'*( X - Q*(Q'*X))
   %            = Q'*X - (Q'*Q) * (Q'*X)
   %            = Q'*X -  (eye) * (Q'*X)
   %            = 0
   % ----------------------------


   % -----------------------------------------
   % note dense blocks are stored in block{} as
   % block * COLUMNS *
   % -----------------------------------------

   for ell=(L-1):-1:(lev+1),
       Q = block{ioff+ell};
       m = size(Q,1);
       mcells = 2^ell;
       isok = (m*mcells == n);
       if (~isok),
         error('OperatorTwoScale_nonwavelet:m=%d,ell=%d,mcells=%d,n=%d', ...
                                               m,   ell,   mcells,   n);
         return;
       end;
         
       % ------------------------------------
       % apply to only the 1st block 
       % since all other blocks are identical
       % ------------------------------------

       for j=1:(mcells/ncells),
         i1 = 1 + (j-1)*m;
         i2 = i1 + m-1;
         QtF = Q(1:m,1:k)'*F( i1:i2, : );
         F(i1:i2,:) = F(i1:i2,:) - Q(1:m,1:k) * QtF(1:k,:);
       end;
   end;
   % -------------------------
   % perform Gram-Schmdit orthogonalization
   % via QR factorization
   % -------------------------
   F = ortho(F);

   % -----------------------------
   % extra the high order terms
   % that have k vanishing moments
   % -----------------------------
   i2 = size(F,2);
   i1 = i2 - k + 1;
   block{ioff+lev} = F(:,i1:i2);
end;

% ------------
% coarse  grid
% ------------
xx = (x - min(x))/(max(x)-min(x));
if (need_shift_x),
  xx = 2*xx - 1;
end;
ndeg = k-1;
F = legendrepoly( ndeg,  x );

% ------------------------------------------------
% apply orthogonalization from all previous levels
% ------------------------------------------------


   for ell=(L-1):-1:0,
       Q = block{ioff + ell};
       m = size(Q,1);
       mcells = 2^ell;
       jsize = n/mcells;

       for j=1:mcells,
         i1 = 1 + (j-1)*m;
         i2 = i1 + m-1;
         QtF = Q(1:m,1:k)'*F( i1:i2, 1:k);
         F(i1:i2,1:k) = F(i1:i2,1:k) - Q(1:m,1:k) * QtF(1:k,1:k);
       end;
   end;

   Q = ortho( F );
   block{ioff + (-1)} = Q;

% ------------------------
% create the sparse matrix
% ------------------------
nzmax = n * k * (L+1);
FMWT = sparse( [], [], [], n,n, nzmax);

% -----------
% Coarse grid
% -----------
irow = 1;
Q = block{ioff+(-1)};
FMWT( irow:(irow+k-1),1:n) = transpose(Q(1:n,1:k));
irow = irow + k;

for lev=0:(L-1),
    Q = block{ioff + lev};
    Qt = transpose( Q );
    ncells = 2^lev;
    isize = n/ncells;
    for i=1:ncells
        i1 = irow;
        i2 = i1 + k-1;

        j1 = 1 + (i-1)*isize;
        j2 = j1 + isize-1;

        FMWT( i1:i2, j1:j2) = Qt(1:k, 1:isize);
        irow = i2 + 1;
    end;
end;



end
