function        [overlap_index] = calc_overlap_index( icell, Lev1, Lev2 );
%
% 
%        [overlap_index] = calc_overlap_index( icell, Lev1, Lev2 );
% 
% compute the overlap_index of cells where icell in Lev1, 
% overlap with cells in Lev2, assume Lev2 >= Lev1
%
% note 0 <= icell <= (2^(Lev1)) - 1
%



overlap_index = [];
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
overlap_index = start_icell:end_icell;

end

