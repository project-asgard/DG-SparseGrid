% simple script to test apply_FMWT_blocks
%
tol = 1e-9;

degree = 3;
level = 2;
max_level = 5;

dof = degree * 2^level;
X = rand(dof, dof);

[FMWT,blocks] = OperatorTwoScale_wavelet2(degree, level);

trans_side = 'LN';
gold = FMWT * X;
assert(test_apply_blocks(gold, X, level, blocks, trans_side, tol));

trans_side = 'RN';
gold = X * FMWT;
assert(test_apply_blocks(gold, X, level, blocks, trans_side, tol));

trans_side = 'LT';
gold = FMWT' * X;
assert(test_apply_blocks(gold, X, level, blocks, trans_side, tol));

trans_side = 'RT';
gold = X * FMWT';
assert(test_apply_blocks(gold, X, level, blocks, trans_side, tol));


function isok = test_apply_blocks(gold, X, level, blocks, ts, tol)
    test = apply_FMWT_blocks(level, blocks, X, ts);
    err = norm(gold - test, 1); 
    rel_err = err / max( norm(gold,1), norm(test,1) );
    if (rel_err >= tol)
       disp(strcat(trans_side,' test failed with relative error ~', num2str(rel_err)));
       isok = false;
    else
       isok = true;
    end
end