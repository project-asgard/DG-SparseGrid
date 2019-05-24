function plot_adapt(pde,pos)

nElements = numel(pde.elementsIDX);
nDims = numel(pde.dimensions);

if nDims == 1
    figure(222);
    subplot(2,3,pos)
    cla
    hold on
    for i=1:nElements
        x = pde.elements.pos_p1(pde.elementsIDX(i))-1;
        y = pde.elements.lev_p1(pde.elementsIDX(i))-1;
        c = pde.elements.node_type(pde.elementsIDX(i));
        if c == 1
            style = 'ob';
        elseif c == 2
            style = 'or';
        end
        offset = 2^(y-1)/2;
        if y > 1
            
            s = 2^(y-1)-1;
            h = 1/(2^(y-1));
            w = 1-h;
            o = h/2;
            plot(x/s*w+o,-y,style);
            
        else
            plot(x+0.5,-y,style);
        end
        xlim([0,1]);
    end
    hold off
end

end