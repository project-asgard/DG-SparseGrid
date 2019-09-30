function tests = md_idx_to_lev_pos_test()
tests = functiontests(localfunctions);
end

function testall(testCase)

max_lev = 8;

%% 1D

gold_lev_vec = [2];
gold_pos_vec = [0];

idx = lev_cell_to_element_index(gold_lev_vec,gold_pos_vec,max_lev);

[lev_vec, pos_vec] = md_idx_to_lev_pos(1, max_lev, idx);

assert(norm(lev_vec-gold_lev_vec) == 0);
assert(norm(pos_vec-gold_pos_vec) == 0);

%% 1D

gold_lev_vec = [0];
gold_pos_vec = [0];

idx = lev_cell_to_element_index(gold_lev_vec,gold_pos_vec,max_lev);

[lev_vec, pos_vec] = md_idx_to_lev_pos(1, max_lev, idx);

assert(norm(lev_vec-gold_lev_vec) == 0);
assert(norm(pos_vec-gold_pos_vec) == 0);

%% 2D

gold_lev_vec = [0,3];
gold_pos_vec = [0,1];

idx = lev_cell_to_element_index(gold_lev_vec,gold_pos_vec,max_lev);

[lev_vec, pos_vec] = md_idx_to_lev_pos(2, max_lev, idx);

assert(norm(lev_vec-gold_lev_vec) == 0);
assert(norm(pos_vec-gold_pos_vec) == 0);

%% 2D

gold_lev_vec = [2,3];
gold_pos_vec = [0,1];

idx = lev_cell_to_element_index(gold_lev_vec,gold_pos_vec,max_lev);

[lev_vec, pos_vec] = md_idx_to_lev_pos(2, max_lev, idx);

assert(norm(lev_vec-gold_lev_vec) == 0);
assert(norm(pos_vec-gold_pos_vec) == 0);

%% 3D

gold_lev_vec = [2,3,4];
gold_pos_vec = [0,1,2];

idx = lev_cell_to_element_index(gold_lev_vec,gold_pos_vec,max_lev);

[lev_vec, pos_vec] = md_idx_to_lev_pos(3, max_lev, idx);

assert(norm(lev_vec-gold_lev_vec) == 0);
assert(norm(pos_vec-gold_pos_vec) == 0);

end