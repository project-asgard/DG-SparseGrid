  % simple test for index_leq_max
% more off;

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
for Levmax=imin:imax,

  result = index_leq_max( ndim, Nmat, Levsum, Levmax );


  m = size(result,1);
  for i=1:m,

     lmax = Nmat{1}( result(i,1));
     lsum = 0;

     for idim=1:ndim,
             iv(idim) = result(i,idim);
             lmax = max(lmax , Nmat{idim}( iv(idim) ));
             lsum = lsum + Nmat{idim}( iv(idim) );
     end;


     isok = (lsum <= Levsum) && (lmax <= Levmax);
     if (~isok),
        % -------------------
        % print debug message
        % -------------------
        msg = sprintf('Levsum=%d,Levmax=%d, ndim=%d,i=%d, ',...
                       Levsum,   Levmax,    ndim,   i);
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
             isok_sum = (Nmat{1}(i1) <= Levsum);
             isok_max = (Nmat{1}(i1) <= Levmax);
             isok = isok_sum && isok_max;  
             if (isok),
                     m = m +1;
             end;
     end;
   elseif (ndim == 2),
     for i1=1:numel( Nmat{1} ),
     for i2=1:numel( Nmat{2} ),
           
           isok_sum =  ((Nmat{1}(i1) + Nmat{2}(i2)) <= Levsum);
           isok_max =  (max(Nmat{1}(i1) , Nmat{2}(i2)) <= Levmax);
           isok = isok_sum && isok_max;
           if (isok),
                   m = m + 1;
           end;
     end;
     end;
   elseif (ndim == 3),
     for i1=1:numel( Nmat{1} ),
     for i2=1:numel( Nmat{2} ),
     for i3=1:numel( Nmat{3} ),
           isok_sum = ((Nmat{1}(i1) + Nmat{2}(i2) + ...
                                      Nmat{3}(i3)) <= Levsum);
           isok_max= (max(Nmat{1}(i1) , ...
                      max(Nmat{2}(i2) , Nmat{3}(i3))) <= Levmax);
           isok = isok_sum && isok_max;
           if (isok),
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
           isok_sum = ((Nmat{1}(i1) + Nmat{2}(i2) + ...
                        Nmat{3}(i3) + Nmat{4}(i4)) <= Levsum);
           isok_max =  (max( max(Nmat{1}(i1) , Nmat{2}(i2)) , ...
                             max(Nmat{3}(i3) , Nmat{4}(i4)) ) <= Levmax);
           isok = isok_sum && isok_max;
           if (isok),
                   m = m + 1;
           end;
     end;
     end;
     end;
     end;
   end;


   isok = (m == size(result,1));
   if (~isok),
     disp(sprintf('Levsum=%d,Levmax=%d, mismatch in m=%d,size(result,1)=%d ', ...
                   Levsum,   Levmax,                m,   size(result,1) )); 
           nerrors = nerrors + 1;
   end;

   if (nerrors == 0),
      disp(sprintf('ndim=%d, Levsum=%d, Levmax=%d is OK', ...
                    ndim,    Levsum,    Levmax));
   else
      error(sprintf('ndim=%d, Levsum=%d, Levmax=%d has %d errors', ...
                     ndim,    Levsum,    Levmax,    nerrors));
      break;
   end;      

end;
end;
end;
                      


% more on;
