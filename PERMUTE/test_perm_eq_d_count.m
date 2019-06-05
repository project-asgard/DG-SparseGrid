% simple test for perm_eq_d_count
nmax = 6;
idim = 2;
for SumEq=0:nmax,
for Lev1=0:nmax,
for Lev2=0:nmax,
   icount = perm_eq_d_count(idim,[Lev1,Lev2], SumEq);
end;
end;
end;


idim = 3;
for SumEq=0:nmax,
for Lev1=0:nmax,
for Lev2=0:nmax,
for Lev3=0:nmax,
   icount = perm_eq_d_count(idim,[Lev1,Lev2,Lev3], SumEq);
end;
end;
end;
end;
