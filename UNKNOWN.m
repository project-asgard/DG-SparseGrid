classdef UNKNOWN < handle
    
    properties
        dimensions = {};
        deg;
        type = 'evolution';
        params = {};
        analytic_solutions = {};
        initial_conditions = {};
        hash_table;
        transform_blocks;
        lo_global; % Start index of unknown in global solution vector
        hi_global; % End   index of unknown in global solution vector
    end
    
    methods
        
        function [ sz ] = size( obj )
            
            num_dims = numel(obj.dimensions);
            
            num_elements = numel(obj.hash_table.elements_idx);
            
            sz = obj.deg^num_dims * num_elements;
            
        end
        
        function fval_initial = get_initial_conditions( obj, opts )
            
            fval_initial...
                = md_eval_function( opts, opts.deg, obj.dimensions,[],...
                                    obj.initial_conditions, obj.hash_table,...
                                    obj.transform_blocks, 0.0 );
            
        end
        
        function unknown = UNKNOWN( opts, dimensions, analytic_solutions, initial_conditions )
            
            unknown.dimensions           = dimensions;
            unknown.deg                  = opts.deg;
            unknown.analytic_solutions   = analytic_solutions;
            unknown.initial_conditions   = initial_conditions;
            
            for i=1:numel(unknown.dimensions)
                unknown.dimensions{i}.lev = opts.lev;
            end
            
            [unknown.hash_table.elements, unknown.hash_table.elements_idx]...
                = hash_table_sparse_nD( unknown.get_lev_vec, opts.max_lev, opts.grid_type );
            
            [~,unknown.transform_blocks] = OperatorTwoScale_wavelet2( opts.deg, opts.max_lev );
            
            compute_dimension_mass_mat( opts, unknown );
            
        end
        
        function lev_vec = get_lev_vec( obj )
            num_dims = numel(obj.dimensions);
            lev_vec  = zeros(num_dims,1);
            for d=1:num_dims
                lev_vec(d,1) = obj.dimensions{d}.lev;
            end
        end
        
        function set_bounds( obj, lo, hi )
            obj.lo_global = lo;
            obj.hi_global = hi;
        end
        
    end
end

