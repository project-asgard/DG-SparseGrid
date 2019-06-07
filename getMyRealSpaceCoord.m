function [coord_D, coord_D_L, coord_D_R] = getMyRealSpaceCoord (pde, opts, hash_table, idx)

%%
% Get the center, left and right coordinates for this element

num_dimensions = numel(pde.dimensions);

for d=1:num_dimensions
    
    xMin = pde.dimensions{d}.domainMin;
    xMax = pde.dimensions{d}.domainMax;
    xRng = xMax-xMin;
    
    %%
    % Get coordinate for this dimension
    
    if opts.use_oldhash
    else
        lev = hash_table.elements.lev_p1(idx,d)-1;
        pos = hash_table.elements.pos_p1(idx,d)-1;
    end
    
    assert(lev <= pde.max_lev);
    assert(lev >= 0);
    assert(pos >= 0);
    
    %%
    % Scaled to [0,1]
    
    if lev > 1
        
        s = 2^(lev-1)-1;
        h = 1/(2^(lev-1));
        w = 1-h;
        o = h/2;
        
        x0 = pos/s*w+o;
%         plot(pos/s*w+o,-lev,style);        
    else
        o = 0.5;
        x0 = pos+0.5;
%         plot(pos+0.5,-lev,style);
    end
    
    %% 
    % Scale to domain
    
    x = x0 * xRng + xMin;
    xL = (x0 - o) * xRng + xMin;
    xR = (x0 + o) * xRng + xMin;
    
    coord_D(d) = x;
    coord_D_L(d) = xL;
    coord_D_R(d) = xR;
    
    assert(xL >= xMin);
    assert(xR <= xMax);
    
end

end