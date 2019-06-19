% simple script to test ConnectnD
%
Dim = 2;
Lev = 6;
nerrors = 0;
for igrid=1:2,
        if (igrid == 1),
                gridType = 'SG';
        else
                gridTye = 'FG';
        end;
%         [HASH, HASHInv] = HashTable(Lev,Dim,gridType);
        lev_vec = zeros(Dim,1) + Lev;
        [HASH, HASHInv] = create_hash_table(lev_vec,gridType);

        Con2D = connect_2D(Lev,HASH,HASHInv,gridType);

        is_sparse_grid =  (strcmp(gridType,'SG'));
        if (is_sparse_grid),
                Levsum = Lev; Levmax = Lev;
        else
                Levsum = Dim*Lev; Levmax = Lev;
        end;
        ConnD = connect_nD(Dim, HASH, HASHInv,Levsum, Levmax);

        % ------------------------
        % compare Con2D and ConnD
        % ------------------------
        isok = (numel(Con2D) == numel(ConnD));
        if (~isok),
          error(sprintf('igrid=%d,Dim=%d,gridType=%s,numel(Con2D),numel(ConnD)', ...
                        igrid,   Dim,   gridType,   numel(Con2D), numel(ConnD)));
          nerrors = nerrors + 1;
        end;

        for i=1:numel(Con2D)
             index_J_2D = Con2D{i};
             index_J_nD = ConnD{i};
             isok = (numel(index_J_2D) == numel(index_J_nD)) && ...
                    all(  sort(index_J_2D) == sort(index_J_nD) );
             if (~isok),
                     nerrors = nerrors + 1;
                     disp('index_J_2D'); sort(index_J_2D);
                     disp('index_J_nD'); sort(index_J_nD);
             end;
        end;
        if (nerrors == 0),
          disp(sprintf('igrid=%d,Lev=%d  is OK',igrid,Lev));
        else
          error(sprintf('igrid=%d,Lev=%d has %d errors',igrid,Lev,nerrors));
        end;
end;
