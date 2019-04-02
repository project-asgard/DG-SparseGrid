function Y = apply_FMWT( kdeg, Lev, FMWT, X, imethod_in, TS_in )
%
% Y = apply_FMWT( k, lev, FMWT, X, [imethod_in], [TS_in] )
% 
% Y = FMWT * X
%
% TS_in =  'LN'  left no-transpose
% TS_in =  'RN'  right no-transpose
% TS_in =  'LT'  left transpose
% TS_in =  'RT'  right transpose
%
idebug = 0;
if (idebug >= 1),
        t1 = cputime();
end;
n = kdeg * 2^Lev;
nrowX = size(X,1);
ncolX = size(X,2);


imethod_default = 3;

imethod = imethod_default;
if (nargin >= 5),
        imethod = imethod_in;
end;

TS_default = 'LN';
TS = TS_default;
if (nargin >= 6),
        TS = TS_in;
end;

if (strcmp(TS,'RN')),
        isLeft = 0; isTrans = 0;
elseif (strcmp(TS,'LT')),
        isLeft = 1; isTrans = 1;
elseif (strcmp(TS,'RT')),
        isLeft = 0; isTrans = 1;
elseif (strcmp(TS,'LN')),
        isLeft = 1; isTrans = 0;
else
   isLeft = 1; isTrans = 0;
end;

isok = (size(FMWT,1) == n) && ...
       (size(FMWT,2) == n);
if (~isok),
  error('apply_FMWT: size(FMWT)=(%d,%d),kdeg=%d,Lev=%d,n=%d', ...
                     size(FMWT,1),size(FMWT,2), kdeg,Lev,n);
  return;
end;

isRight = ~isLeft;

% -------------------------
% determine shape of output
% -------------------------
   isok = (isLeft  && (nrowX == n)) || ...
          (isRight && (ncolX == n));
   if (~isok),
      error('apply_FMWT: isLeft=%d, nrowX=%d,ncolX=%d,n=%d', ...
                         isLeft,    nrowX,   ncolX,   n);
      return;
   end;


   if (isLeft),
     nrowY = n; ncolY = ncolX;
   else
     nrowY = nrowX; ncolY = n;
   end;
   Y = zeros(nrowY, ncolY);


if (imethod == 1),
      % ----------------------------
      % just perform matrix multiply
      % ----------------------------
      if (isLeft),
        if (isTrans),
         Y = FMWT'*X;
        else
         Y = FMWT * X;
        end;
      else
        if (isTrans),
           Y = X * FMWT';
        else
           Y = X * FMWT;
        end;
      end;
elseif (imethod == 2),
      % -------------------------------------------- 
      % take advantage of special sparsity structure
      % -------------------------------------------- 

      % -------------------------------
      % special case for coarsest level
      % -------------------------------
      ip = 1;
      ipend = 2*kdeg;
      col1 = 1;
      col2 = n;

      Fmat = full(FMWT(ip:ipend, col1:col2));
      if (isLeft),
         if (isTrans),
           Y(col1:col2,1:ncolX) = Fmat' * X( ip:ipend,1:ncolX);
         else
           Y(ip:ipend,1:ncolX) = Fmat * X( col1:col2,1:ncolX);
         end;
      else
         if (isTrans),
           Y(1:nrowX,ip:ipend) = X(1:nrowX, col1:col2) * Fmat';
         else
           Y(1:nrowX,col1:col2) = X(1:nrowX, ip:ipend) * Fmat;
         end;
      end;

      ip = 2*kdeg + 1;
      for lev=1:(Lev-1),
          ncells = 2^lev;
          isize = n/ncells;
          for icell=1:ncells,
             ipend = ip + kdeg-1;
                  col1 = 1 + (icell-1)*isize;
                  col2 = col1 + isize-1;
		  Fmat = full( FMWT(ip:ipend, col1:col2) );
                  if (isLeft),
                     if (isTrans),
                       Y(col1:col2,1:ncolX) = Y(col1:col2,1:ncolX) + ...
                                              Fmat' * X( ip:ipend,1:ncolX);
                     else
                       Y(ip:ipend,1:ncolX) = Y(ip:ipend,1:ncolX) + ...
                                             Fmat * X( col1:col2,1:ncolX);
                     end;
                  else
                     if (isTrans),
                       Y(1:nrowX,ip:ipend) = Y(1:nrowX,ip:ipend) + ...
                                             X(1:nrowX, col1:col2) * Fmat'; 
                     else
                       Y(1:nrowX,col1:col2) = Y(1:nrowX,col1:col2) + ...
                                              X(1:nrowX, ip:ipend) * Fmat; 
                     end;
                  end;

             ip = ipend + 1;
          end;
      end;
elseif (imethod == 3),

      % -------------------------------------------- 
      % take advantage of special sparsity structure
      % and identical sub-matrices
      % -------------------------------------------- 

      % -------------------------------
      % special case for coarsest level
      % -------------------------------
      ip = 1;
      ipend = 2*kdeg;
      col1 = 1;
      col2 = n;
      Fmat = full( FMWT(ip:ipend,col1:col2) );
      if (isLeft),
         if (isTrans),
            Y(col1:col2,1:ncolX) = Fmat' * X( ip:ipend,1:ncolX);
         else
            Y(ip:ipend,1:ncolX) = Fmat * X( col1:col2,1:ncolX);
         end;
      else
         if (isTrans),
            Y(1:nrowX,ip:ipend) = X(1:nrowX, col1:col2) * Fmat'; 
         else
            Y(1:nrowX,col1:col2) = X(1:nrowX, ip:ipend) * Fmat;
         end;
      end;

      ip = 2*kdeg + 1;
      for lev=1:(Lev-1),
          ncells = 2^lev;
          isize = n/ncells;

          % --------------------------------------
          % take advantage of identical sub-blocks
          % to make fewer calls to batched GEMM
          % --------------------------------------
          icell = 1;
          ipend = ip + kdeg-1;
          col1 = 1 + (icell-1)*isize;
          col2 = col1 + isize-1;

          Fmat = full(FMWT(ip:ipend,col1:col2));
          ip2 = ip + (ncells * isize)-1;

          ncol = (n*ncolX/isize);
          nrows = (ncells*kdeg);
          ipend = ip + (ncells*kdeg)-1;

          if (isLeft),
            if (isTrans),
              Y(1:nrowY,1:ncolY) = Y(1:nrowY,1:ncolY) + ...
                  reshape( ...
                    Fmat(1:kdeg,1:isize)' * ...
                        reshape( X(ip:ipend,1:ncolX), kdeg, ((ipend-ip+1)*ncolX/kdeg)), ...
                    nrowY, ncolY );
            else
             Y(ip:ipend,1:ncolX) = ...
               reshape(  Fmat(1:kdeg,1:isize)*reshape(X(1:n,1:ncolX), isize, ncol), ...
                     nrows, (kdeg*ncol)/nrows  );
            end;
          else
            if (isTrans),
              for icell=1:ncells,
                 j1 = ip + (icell-1)*kdeg;
                 j2 = j1 + kdeg-1;

                 jx1 = 1 + (icell-1)*isize;
                 jx2 = jx1 + isize-1;
                 Y(1:nrowX, j1:j2) = ...
                     X(1:nrowX, jx1:jx2) * Fmat(1:kdeg,1:isize)';
               end; % end for
            else
              for icell=1:ncells,
                  j1 = 1 + (icell-1)*isize;
                  j2 = j1 + isize-1;

                  jx1 = ip + (icell-1)*kdeg;
                  jx2 = jx1 + kdeg-1;
                  Y(1:nrowX, j1:j2) = Y(1:nrowX, j1:j2) + ...
                       X(1:nrowX,  jx1:jx2 ) * Fmat( 1:kdeg, 1:isize );
               end; % end for
            end;
          end;


          ip = ip + (ncells * kdeg);

      end;


else
    error('apply_FMWT: invalid imethod=%d',imethod);
    return;
end;

if (idebug >= 1),
   t2 = cputime();
   disp(sprintf('apply_FMWT: n=%d, nrowX=%d, ncolX=%d, imethod=%d, time=%g', ...
        n,  nrowX, ncolX, imethod, (t2-t1) ));
end;

end

   

