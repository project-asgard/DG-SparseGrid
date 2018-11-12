% simple test for index_max

Nmat{1} = 2:3;
Nmat{2} = 0:4;
Nmat{3} = 0:3;
Nmat{4} = 1:5;

for i=1:numel(Nmat),
        minval(i) = min( Nmat{i} );
        maxval(i) = max( Nmat{i} );
end;

nerrors = 0;

% ----------------
% N dimension case
% ----------------
for ndim = 1:4,
imin = min( minval(1:ndim) );
imax = sum( maxval(1:ndim) );
for Levmax=imin:imax,

  result = index_max( ndim, Nmat, Levmax );


  m = size(result,1);
  for i=1:m,

     lmax = Nmat{1}( result(i,1));
     for idim=1:ndim,
             iv(idim) = result(i,idim);
             lmax = max(lmax , Nmat{idim}( iv(idim) ));
     end;


     isok = (lmax <= Levmax);
     if (~isok),
        % -------------------
        % print debug message
        % -------------------
        msg = sprintf('Levmax=%d, ndim=%d,i=%d, ',...
                       Levmax,ndim,i);
        for idim=1:ndim,
          msg = strcat(msg,sprintf('Nmat{%d}(%d)=%d, ', ...
                        idim, iv(idim), Nmat{idim}(iv(idim)) ));
        end;
        disp( msg );
        nerrors = nerrors + 1;
     end;
   end;

   m = 0;
   if (ndim == 1),
     for i1=1:numel( Nmat{1} ),
             if (Nmat{1}(i1) <= Levmax),
                     m = m +1;
             end;
     end;
   elseif (ndim == 2),
     for i1=1:numel( Nmat{1} ),
     for i2=1:numel( Nmat{2} ),
           if (max(Nmat{1}(i1) , Nmat{2}(i2)) <= Levmax),
                   m = m + 1;
           end;
     end;
     end;
   elseif (ndim == 3),
     for i1=1:numel( Nmat{1} ),
     for i2=1:numel( Nmat{2} ),
     for i3=1:numel( Nmat{3} ),
           if (max(Nmat{1}(i1) , max(Nmat{2}(i2) , Nmat{3}(i3))) <= Levmax),
                   m = m + 1;
           end;
     end;
     end;
     end;
   elseif (ndim == 4),
     for i1=1:numel( Nmat{1} ),
     for i2=1:numel( Nmat{2} ),
     for i3=1:numel( Nmat{3} ),
     for i4=1:numel( Nmat{4} ),
           if (max( max(Nmat{1}(i1) , Nmat{2}(i2)) , ...
                    max(Nmat{3}(i3) , Nmat{4}(i4)) ) <= Levmax),
                   m = m + 1;
           end;
     end;
     end;
     end;
     end;
   end;


   isok = (m == size(result,1));
   if (~isok),
           disp(sprintf('mismatch in m=%d,size(result,1)=%d ', ...
                                     m,   size(result,1) )); 
           nerrors = nerrors + 1;
   end;

   if (nerrors == 0),
      disp(sprintf('ndim=%d, Levmax=%d is OK', ndim,Levmax));
   else
      error(sprintf('ndim=%d, Levmax=%d has %d errors', ...
                    ndim,    Levmax, nerrors));
      break;
   end;      

end;
end;
                      


