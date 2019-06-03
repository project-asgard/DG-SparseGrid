function plot_adapt(pde,pos)

num_elements = numel(pde.elementsIDX);
num_dims = numel(pde.dimensions);
xMin = pde.dimensions{1}.domainMin;
xMax = pde.dimensions{1}.domainMax;


%%
% Get grid point coordinates

coordinates = zeros (num_elements, num_dims);
levels      = zeros (num_elements, 1);
style       = strings (num_elements, 1);

for i=1:num_elements
    
    idx = pde.elementsIDX(i);
    c = pde.elements.node_type(idx);
    
    if c == 1
        style(i) = 'ob';
    elseif c == 2
        style(i) = 'or';
    end
    
    coordinates(i,:) = getMyRealSpaceCoord(pde,idx);
    
    levels(i) = sum (pde.elements.lev_p1(idx,:)-1,'all');

end

figure(222);
subplot(2,3,pos)
cla
hold on

if num_dims == 1

    for i=1:num_elements
        
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

elseif num_dims == 2

    for i=1:num_elements
        plot3 (coordinates(i,1), coordinates(i,2), levels(i), style(i));
        hold on
    end

end

hold off

end