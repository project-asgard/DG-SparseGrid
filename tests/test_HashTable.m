% more off;
gridType_table{1} = 'SG';
gridType_table{2} = 'FG';
ncase = numel(gridType_table);

nerrors = 0;
Lev = 4;
Dim = 2;
for icase=1:ncase,
  gridType = gridType_table{icase}; 
  time_2D = tic();
%   forwardHash2D,inverseHash2D] = HashTable2D(Lev,Dim,gridType);
  [forwardHash2D,inverseHash2D] = hash_table_2D(Lev,Dim,gridType);
  elapsed_time_2D = toc( time_2D );

  time_nD = tic();
%   [forwardHash, inverseHash] = HashTable(Lev,Dim,gridType);
  lev_vec = zeros(Dim,1)+Lev;
  [forwardHash, inverseHash] = hash_table_nD(lev_vec,gridType);
  elapsed_time_nD = toc( time_nD );

  disp(sprintf('Lev=%d,Dim=%d, time for HashTable2D=%g, time for HashTable=%g',...
                Lev,   Dim,    elapsed_time_2D, elapsed_time_nD ));
  
  disp(sprintf('numel(inverseHash)=%g', numel(inverseHash)));

  isok = numel(fieldnames( forwardHash2D)) == numel(fieldnames(forwardHash)) && ...
         numel( inverseHash2D) == numel(inverseHash);
  if (~isok),
     disp(sprintf('numel(fieldnames(forwardHash2D)) = %d numel(fieldnames(forwardHash)) = %d', ...
                   numel(fieldnames(forwardHash2D)),     numel(fieldnames(forwardHash)) ));
  
    
     disp(sprintf('numel(inverseHash2D) = %d numel(inverseHash) = %d', ...
                   numel(inverseHash2D),     numel(inverseHash) ));
     nerrors = nerrors + 1;
  end;
  
  if (~exist('hash_format')),
      %     hash_format = 'i%04.4d';
      set_hash_format
  end;
  
  % -------------------------------------------------
  % check key in forwardHash2D is in other Hash table
  % -------------------------------------------------
  
  nelem = numel(inverseHash2D);
  is_seen = zeros( nelem,1);
  
  for icell=1:numel(inverseHash2D),
          ival = inverseHash2D{icell};
          key_ival = ival(1:4);
          key_string = sprintf(hash_format,key_ival);
          isok = isfield( forwardHash, key_string);
          if (~isok),
                  disp(sprintf('icell=%d,key_string=%s', ...
                                icell,   key_string));
                  nerrors = nerrors + 1;
          end;
          icount = getfield(forwardHash,key_string);
          is_seen(icount) = is_seen(icount) + 1;
  end;
  isok = all( is_seen(1:nelem) == ones(nelem,1) );
  if (~isok),
          idx = find( is_seen ~= 1 );
          disp('idx')
          idx
          disp('is_seen(idx)')
          is_seen(idx)

          nerrors  = nerrors + 1;
  end;
  
  
  % -------------------------------------------------
  % check key in forwardHash is in other Hash table
  % -------------------------------------------------
  is_seen = zeros( nelem,1);
  
  for icell=1:numel(inverseHash),
          ival = inverseHash{icell};
          key_ival = ival(1:4);
          key_string = sprintf(hash_format,key_ival);
          isok = isfield( forwardHash2D, key_string);
          if (~isok),
                  disp(sprintf('icell=%d,key_string=%s', ...
                                icell,   key_string));

                  nerrors = nerrors + 1;
          end;
          icount = getfield(forwardHash2D,key_string);
          is_seen(icount) = is_seen(icount) + 1;
  end;
  isok = all( is_seen(1:nelem) == ones(nelem,1) );
  if (~isok),
          idx = find( is_seen ~= 1 );
          disp('idx')
          idx
          disp('is_seen(idx)')
          is_seen(idx)

          nerrors = nerrors + 1;
  end;
  
end;  

if (nerrors == 0),
        disp('ALL OK ');
end;

% -----------------------------
% test performance of HashTable, Dim=4
% -----------------------------
gridType = 'SG';
Dim = 4;
maxLev = 8;
elapsed_time = zeros(maxLev,1);
sizes = zeros(maxLev,1);
for Lev=1:maxLev,
   t1 = tic();
%    [forwardHash,inverseHash] = HashTable(Lev,Dim,gridType);
   lev_vec = zeros(Dim,1)+Lev;
   [forwardHash,inverseHash] = hash_table_nD(lev_vec,gridType);
   elapsed_time(Lev) = toc( t1 );

   sizes(Lev) = numel(inverseHash);
end;

clf;
figure(1);
subplot(2,1,1); plot( 1:maxLev, elapsed_time(1:maxLev) ); 
title(sprintf('elapsed time, Dim=%d',Dim));
subplot(2,1,2); plot( 1:maxLev, sizes(1:maxLev) );
title(sprintf('size of inverseHash, Dim=%d',Dim));
% print -djpg test_HashTable_Dim=4.jpg
filename = sprintf('test_HashTable_Dim=%d.jpg', Dim);
print(filename,'-djpeg');




% -----------------------------
% test performance of HashTable, Dim=6
% -----------------------------
gridType = 'SG';
Dim = 6;
maxLev = 8;
elapsed_time = zeros(maxLev,1);
sizes = zeros(maxLev,1);
for Lev=1:maxLev,
   t1 = tic();
%    [forwardHash,inverseHash] = HashTable(Lev,Dim,gridType);
    lev_vec = zeros(Dim,1)+Lev;
    [forwardHash,inverseHash] = hash_table_nD(lev_vec,gridType);

   elapsed_time(Lev) = toc( t1 );

   sizes(Lev) = numel(inverseHash);
end;

clf;
figure(2);
subplot(2,1,1); plot( 1:maxLev, elapsed_time(1:maxLev) ); 
title(sprintf('elapsed time, Dim=%d',Dim));
subplot(2,1,2); plot( 1:maxLev, sizes(1:maxLev) );
title(sprintf('size of inverseHash, Dim=%d',Dim));
% print -djpg test_HashTable_Dim=6.jpg
filename = sprintf('test_HashTable_Dim=%d.jpg', Dim);
print(filename,'-djpeg');


% more on;
