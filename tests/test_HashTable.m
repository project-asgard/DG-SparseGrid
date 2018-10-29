gridType_table{1} = 'SG';
gridType_table{2} = 'FG';
ncase = numel(gridType_table);

Lev = 4;
Dim = 2;
for icase=1:ncase,
  gridType = gridType_table{icase}; 
  [forwardHash2D,inverseHash2D] = HashTable2D(Lev,Dim,gridType);
  [forwardHash, inverseHash] = HashTable(Lev,Dim,gridType);
  
  isok = numfields( forwardHash2D) == numfields(forwardHash) && ...
         numel( inverseHash2D) == numel(inverseHash);
  if (~isok),
     disp(sprintf('numfields(forwardHash2D) = %d numfields(forwardHash) = %d', ...
                   numfields(forwardHash2D),     numfields(forwardHash) ));
  
    
     disp(sprintf('numel(inverseHash2D) = %d numel(inverseHash) = %d', ...
                   numel(inverseHash2D),     numel(inverseHash) ));
  end;
  
  if (~exist('hash_format')),
    hash_format = 'i%04.4d_';
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
  end;
  
end;  
