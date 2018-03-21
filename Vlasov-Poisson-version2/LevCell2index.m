% for given Lev and Cell, determine Index
function index=LevCell2index(Lev,Cell)

index=2.^(Lev-1)+Cell+1;

ix = find(Lev ==0);
index(ix)=1;

end