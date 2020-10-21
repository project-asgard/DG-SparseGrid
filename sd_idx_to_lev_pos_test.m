function tests = sd_idx_to_lev_pos_test
tests = functiontests(localfunctions);
end


function sd_test(testCase)

[lev,pos] = sd_idx_to_lev_pos (1);

gold_lev = 0;           
gold_pos = 0;

assert( lev == gold_lev );
assert( pos == gold_pos );

[lev,pos] = sd_idx_to_lev_pos (2);

gold_lev = 1;           
gold_pos = 0;

assert( lev == gold_lev );
assert( pos == gold_pos );

[lev,pos] = sd_idx_to_lev_pos (5);

gold_lev = 3;           
gold_pos = 0;

assert( lev == gold_lev );
assert( pos == gold_pos );

[lev,pos] = sd_idx_to_lev_pos (6);

gold_lev = 3;           
gold_pos = 1;

assert( lev == gold_lev );
assert( pos == gold_pos );


[lev,pos] = sd_idx_to_lev_pos (7);

gold_lev = 3;           
gold_pos = 2;

assert( lev == gold_lev );
assert( pos == gold_pos );

[lev,pos] = sd_idx_to_lev_pos (8);

gold_lev = 3;           
gold_pos = 3;

assert( lev == gold_lev );
assert( pos == gold_pos );

end


function sd2_test(testCase)

%%

gold_lev = 0;           
gold_pos = 0;

idx = lev_cell_to_1D_index(gold_lev,gold_pos);

[lev,pos] = sd_idx_to_lev_pos(idx);

assert( lev == gold_lev );
assert( pos == gold_pos );

%%

gold_lev = 1;           
gold_pos = 0;

idx = lev_cell_to_1D_index(gold_lev,gold_pos);

[lev,pos] = sd_idx_to_lev_pos(idx);

assert( lev == gold_lev );
assert( pos == gold_pos );

%%

gold_lev = 3;           
gold_pos = 0;

idx = lev_cell_to_1D_index(gold_lev,gold_pos);

[lev,pos] = sd_idx_to_lev_pos(idx);

assert( lev == gold_lev );
assert( pos == gold_pos );

%%

gold_lev = 5;           
gold_pos = 5;

idx = lev_cell_to_1D_index(gold_lev,gold_pos);

[lev,pos] = sd_idx_to_lev_pos(idx);

assert( lev == gold_lev );
assert( pos == gold_pos );

end
