function [Y] = get_moments_for_plot(pde)
%Gets the moment values used in plotting.  Outputs a matrix that has
%num_moments columns and each column is the moment in time.


num_moments = numel(pde.moments);
Y = [];

for m=1:num_moments
    moment = pde.moments{m};
   
    Y = [Y moment.moment_fval_integral'];
end

end
