function [pde, moment_values, moment_analytic_values] = calculate_mass(pde,opts,fval,hash_table,t)

num_moments = numel(pde.moments);
moment_values = zeros(num_moments,1);
moment_analytic_values = [];

%Find analytic integral
if ~isempty(pde.solutions)
    %Get solution
    mass_bool = 1; %Need actual solution.
    soln = md_eval_function(opts,opts.deg,pde.dimensions,pde.params,pde.solutions,hash_table,pde.transform_blocks,t,mass_bool);     
    moment_analytic_values = zeros(num_moments,1);
end

for m=1:num_moments
    moment = pde.moments{m};
    %Calculate vector
    moment = moment.createMomentVector(opts,hash_table);
    
    %Calculate actual moment through dot product
    moment_values(m) = moment.vector'*fval;
    
    %Store with previous moments
    %moment.moment_fval_integral = [moment.moment_fval_integral moment_values(m)];
    
    if ~isempty(pde.solutions)
        %Calculate actual moment through dot product
        moment_analytic_values(m) = moment.vector'*soln;
        
        %Store with previous moments
        %moment.moment_analytic_integral = [moment.moment_analytic_integral ...
        %                    moment_analytic_values(m)];
    end
    
    pde.moments{m} = moment;
end

end