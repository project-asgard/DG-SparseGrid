function coord = get_realspace_coords(pde,nodes)

%%
% Yes, this is crap. 

nDims = numel(pde.dimensions);

coord = [];

if nDims <= 3
    
    if nDims ==1
        [xx1] = ndgrid(nodes{1});
        coord = {xx1};
    end
    if nDims==2
        [xx1,xx2] = ndgrid(nodes{2},nodes{1});
        coord = {xx2,xx1};
    end
    if nDims==3
        [xx1,xx2,xx3] = ndgrid(nodes{3},nodes{2},nodes{1});
        coord = {xx3,xx2,xx1};
    end
    if nDims==6
        [xx1,xx2,xx3,xx4,xx5,xx6] = ndgrid(nodes{6},nodes{5},nodes{4},nodes{3},nodes{2},nodes{1});
        coord = {xx6,xx5,xx4,xx3,xx2,xx1};
    end
    
end

end