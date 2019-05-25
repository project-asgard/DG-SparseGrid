function plot_adapt(pde,pos)

nElements = numel(pde.elementsIDX);
nDims = numel(pde.dimensions);
xMin = pde.dimensions{1}.domainMin;
xMax = pde.dimensions{1}.domainMax;

if nDims == 1
    figure(222);
    subplot(2,3,pos)
    cla
    hold on
    for i=1:nElements
        
        idx = pde.elementsIDX(i);
        y = pde.elements.lev_p1(idx,1)-1;
        c = pde.elements.node_type(idx);
        if c == 1
            style = 'ob';
        elseif c == 2
            style = 'or';
        end
        
        coord_D = getMyRealSpaceCoord(pde,idx);
        plot(coord_D(1),-y,style);
        
        xlim([xMin,xMax]);
    end
    hold off
end

end