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
        
        function [ sz ] = size( obj, use_oldhash )
            
            if( nargin < 2 )
                use_oldhash = false;
            end
            
            num_dims = numel(obj.dimensions);
            
            if( use_oldhash )
                num_elements = numel(obj.hash_table);
            else
                num_elements = numel(obj.hash_table.elements_idx);
            end
            
            sz = obj.deg^num_dims * num_elements;
            
        end
        
        function set_initial_conditions( obj, opts )
            
            obj.fval = md_eval_function( opts, opts.deg, obj.dimensions,[],...
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
        
        function lev_vec = get_lev_vec(obj)
            num_dims = numel(obj.dimensions);
            lev_vec  = zeros(num_dims,1);
            for d=1:num_dims
                lev_vec(d,1) = obj.dimensions{d}.lev;
            end
        end
        
    end
end

