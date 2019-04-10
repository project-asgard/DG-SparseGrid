% simple test for index_leq

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
for Levsum=imin:imax,

  result = index_leq( ndim, Nmat, Levsum );


  m = size(result,1);
  for i=1:m,
     isum = 0;
     for idim=1:ndim,
             iv(idim) = result(i,idim);
             isum = isum + Nmat{idim}( iv(idim) );
     end;


     isok = (isum <= Levsum);
     if (~isok),
        % -------------------
        % print debug message
        % -------------------
        msg = sprintf('ndim=%d,i=%d,',ndim,i);
        for idim=1:ndim,
          msg = strcat(msg,sprintf('Nmat{%d}(%d)=%d,', ...
                        idim, iv(idim), Nmat{idim}(iv(idim)) ));
        end;
        disp( msg );
        nerrors = nerrors + 1;
     end;
   end;

   m = 0;
   if (ndim == 1),
     for i1=1:numel( Nmat{1} ),
             if (Nmat{1}(i1) <= Levsum),
                     m = m +1;
             end;
     end;
   elseif (ndim == 2),
     for i1=1:numel( Nmat{1} ),
     for i2=1:numel( Nmat{2} ),
           if (Nmat{1}(i1) + Nmat{2}(i2) <= Levsum),
                   m = m + 1;
           end;
     end;
     end;
   elseif (ndim == 3),
     for i1=1:numel( Nmat{1} ),
     for i2=1:numel( Nmat{2} ),
     for i3=1:numel( Nmat{3} ),
           if (Nmat{1}(i1) + Nmat{2}(i2) + ...
               Nmat{3}(i3) <= Levsum),
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
           if (Nmat{1}(i1) + Nmat{2}(i2) + ...
               Nmat{3}(i3) + Nmat{4}(i4) <= Levsum),
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
      disp(sprintf('ndim=%d, Levsum=%d is OK', ndim,Levsum));
   else
      error(sprintf('ndim=%d, Levsum=%d has %d errors', ...
                    ndim,    Levsum, nerrors));
      break;
   end;      

end;
end;
                      


