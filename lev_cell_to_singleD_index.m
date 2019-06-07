function index = lev_cell_to_singleD_index(lev,cell)

%%
% Map lev,cell to index within a single dim

index=2.^(lev-1)+cell+1;

ix = find(lev == 0);
index(ix)=1;

end