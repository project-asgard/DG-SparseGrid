function [lev,pos] = sd_idx_to_lev_pos (idx)

debug = 0;

if debug; disp(idx); end

assert(idx >= 1);
assert(fix(idx)==idx);

if idx == 1
    
    lev = 0;
    pos = 0;
    
else
    
    assert(idx-1<=1e19);
    
    lev = floor(log2(double(idx-1))+1);
    pos = double((idx-1)-2^(lev-1));
    
end

assert(lev>=0);
assert(pos>=0);

end
