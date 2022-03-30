classdef PDE_SOLUTIONS < handle
    %UNTITLED4 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        dimensions = {};
        num_funcs
        analytic_solutions = {};
        initial_conditions = {};
        fval;
        F_fval; % F(fval)
        transform_blocks;
        params = {};
        %hash_table?
    end
    
    methods
        
        function set_initial_conditions( obj, opts, hash_table )
            
            for i = 1 : obj.num_funcs
                
                obj.fval(:,i)...
                    = md_eval_function( opts, opts.deg, obj.dimensions,[],...
                                        obj.initial_conditions(i), hash_table, obj.transform_blocks, 0.0 );
                
            end
            
        end
        
        function initialize_auxiliary( obj )
            
            obj.F_fval = zeros( size( obj.fval ) );
            
        end
        
        function pde_solutions = PDE_SOLUTIONS( opts, dimensions, analytic_solutions, initial_conditions )
            
            [~,pde_solutions.transform_blocks] = OperatorTwoScale_wavelet2(opts.deg, opts.max_lev);
            pde_solutions.num_funcs            = numel(initial_conditions);
            pde_solutions.dimensions           = dimensions;
            pde_solutions.analytic_solutions   = analytic_solutions;
            pde_solutions.initial_conditions   = initial_conditions;
            
            pde_solutions = compute_dimension_mass_mat( opts, pde_solutions );
            
        end
        
    end
end

