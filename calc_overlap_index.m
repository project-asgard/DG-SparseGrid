function        [start_icell, end_icell] = calc_overlap_index( n1, i1, Lev, Deg );
%
% 
%        [start_icell, end_icell] = calc_overlap_index( n1,i1,Lev,Deg
% 
% compute the overlap_index of cells where icell in Lev1, 
% overlap with cells in Lev2, assume Lev2 >= Lev1
%
% note 0 <= icell <= (2^(Lev1)) - 1
%
icell = i1;
Lev1 = n1;
Lev2 = Lev;



start_icell = 0;
end_icell = 0;
% ----------------------------------------------
% note that each cell is subdivided by 2^(Lev2-Lev1) times
% ----------------------------------------------
isok_icell = (0 <= icell) && (icell <= (2^Lev1)-1);
if (~isok_icell),
  error(sprintf('calc_overlap_index: invalid icell=%d, Lev1=%d',  ...
                                             icell,    Lev1 ));
  return;
end;

isok_Lev2 = (0 <= Lev1) && (Lev1 <= Lev2);
if (~isok_Lev2),
  error(sprintf('calc_overlap_index: invalid Lev2=%d, Lev1=%d', ...
                                             Lev2,    Lev1 ));
  return;
end;

isize = 2^(Lev2-Lev1);
start_icell = icell * 2^(Lev2-Lev1);
end_icell = start_icell + isize - 1;

csize = 2*(Deg-1);
start_icell = start_icell * csize + 1;
end_icell = end_icell * csize + 1 + (csize-1);
end_icell  = min(end_icell, 2^(Lev-1) * csize );

end

