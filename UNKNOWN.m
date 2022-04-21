classdef UNKNOWN < handle
    
    properties
        dimensions = {};
        num_funcs;
        deg;
        params = {};
        analytic_solutions = {};
        initial_conditions = {};
        fval;
        hash_table;
        transform_blocks;
    end
    
    methods
        
        function set_initial_conditions( obj, opts )
            
            for i = 1 : obj.num_funcs
                
                obj.fval(:,i)...
                    = md_eval_function( opts, opts.deg, obj.dimensions,[],...
                                        obj.initial_conditions(i), obj.hash_table, obj.transform_blocks, 0.0 );
                
            end
            
        end
        
        function unknown = UNKNOWN( opts, dimensions, num_funcs, analytic_solutions, initial_conditions )
            
            unknown.dimensions           = dimensions;
            unknown.num_funcs            = num_funcs;
            unknown.deg                = opts.deg;
            unknown.analytic_solutions   = analytic_solutions;
            unknown.initial_conditions   = initial_conditions;
            
            [unknown.hash_table.elements, unknown.hash_table.elements_idx]...
                = hash_table_sparse_nD( unknown.get_lev_vec, opts.max_lev, opts.grid_type );
            
            [~,unknown.transform_blocks] = OperatorTwoScale_wavelet2( opts.deg, opts.max_lev );
            
            compute_dimension_mass_mat( opts, unknown );
            
        end
        
        function lev_vec = get_lev_vec(obj)
            num_dims = numel(obj.dimensions);
            lev_vec  = zeros(num_dims,1);
            for d=1:num_dims
                lev_vec(d,1) = obj.dimensions{d}.lev;
            end
        end
        
    end
end

