

disp('simple test of perm_leq_max ');
nerrors = 0;

idim = 6;
n_leq = 6;
n_max = 4;

disp(sprintf('idim=%d, n_leq=%d, n_max=%d', ...
              idim, n_leq, n_max ));


for order_by_n=1:-1:0,
   disp(sprintf('order_by_n=%d', order_by_n));

   result = perm_leq_max( idim, n_leq, n_max, order_by_n );
  
   disp(sprintf('size(result,1) = %g ', size(result,1)));
   
  

   isok = all( max( result, [], 2) <= n_max ) && ...
          all( sum( result, 2) <= n_leq );
   if (~isok),
           idx = find( (max(result,[],2) > n_max) ||  ...
                       (sum( result > n_leq))  );
           for i=1:numel(idx),
                   j = idx(i);
                   disp(sprintf('j=%d,result(j,:) ',j));
                   result(j,1:idim)
           end;
           nerrors = nerrors + 1;
   end;
end;

if (nerrors == 0),
        disp('ALL OK');
end;
