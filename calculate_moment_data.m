function [pde] = calculate_moment_data(pde,opts)

num_moments = numel(pde.moments);

for m=1:num_moments
    %Calculate vector
    pde.moments{m} = pde.moments{m}.createFlist(pde,opts);
end

end